"""Small fixed geometry/network scores. No training or target geometry input."""
import csv
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from scipy.stats import rankdata, spearmanr

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
from pilot import read_star, metric, sha

PROTOCOL_SHA='f41cf221d1d3e798fb6708ef3c326ea6b505eb829a718d16e68b4d94a3d170d0'
PRIMARY='sidechain_amide_distance'
SCORES=['geometry_locality',PRIMARY,'sidechain_ca_distance','all_heavy_distance','contact_diffusion']


def minimum_distance(a,b):
    a,b=np.asarray(a,float),np.asarray(b,float)
    return float(np.linalg.norm(a[:,None,:]-b[None,:,:],axis=-1).min())


def diffuse(ca,source,restart=.5):
    d=np.linalg.norm(ca[:,None,:]-ca[None,:,:],axis=-1)
    a=((d<=8)&(d>0)).astype(float)
    for i in np.flatnonzero(a.sum(axis=1)==0):a[i,i]=1.
    p=a/a.sum(axis=1,keepdims=True)
    seed=np.zeros(len(ca));seed[source]=restart
    result=np.linalg.solve(np.eye(len(ca))-(1-restart)*p.T,seed)
    assert np.all(result>=-1e-12) and np.isclose(result.sum(),1)
    return result


def source_scores(path):
    # Parse only the pre-mutation/double-mutant structure.
    chain=next(MMCIFParser(QUIET=True).get_structure('source',str(path))[0].get_chains())
    residues={r.id[1]:r for r in chain if r.id[0]==' ' and 'CA' in r}
    assert residues[97].resname=='ILE','Source must contain I97, not target V97'
    side=np.array([a.coord for a in residues[97] if a.element!='H' and a.name not in ('N','CA','C','O','OXT')])
    assert len(side)==4, 'Expected I97 sidechain CB, CG1, CG2, CD1'
    ids=sorted(residues)
    ca=np.array([residues[i]['CA'].coord for i in ids],float)
    walks=diffuse(ca,ids.index(97))
    result={}
    for j,pos in enumerate(ids):
        r=residues[pos]
        if 'N' not in r:continue
        heavy=np.array([a.coord for a in r if a.element!='H'])
        result[pos]={
            'geometry_locality':-float(np.linalg.norm(ca[j]-ca[ids.index(97)])),
            PRIMARY:-minimum_distance(side,[r['N'].coord]),
            'sidechain_ca_distance':-minimum_distance(side,[r['CA'].coord]),
            'all_heavy_distance':-minimum_distance(side,heavy),
            'contact_diffusion':float(walks[j])}
    return result


def rho(x,y):
    a,b=rankdata(x),rankdata(y)
    return float(np.corrcoef(a,b)[0,1]) if np.ptp(a)>0 and np.ptp(b)>0 else None


def bootstrap(rows):
    # Rank data are paired within each resample. Sequence blocks retain some local dependence.
    rng=np.random.default_rng(20261110)
    n=len(rows)
    a=np.array([r[PRIMARY] for r in rows]);b=np.array([r['geometry_locality'] for r in rows]);y=np.array([r['csp'] for r in rows])
    delta=[]
    for _ in range(2000):
        starts=rng.integers(n,size=math.ceil(n/10))
        ix=np.concatenate([(start+np.arange(10))%n for start in starts])[:n]
        ra,rb=rho(a[ix],y[ix]),rho(b[ix],y[ix])
        if ra is not None and rb is not None:delta.append(ra-rb)
    return {'resamples':len(delta),'circular_block_length':10,'seed':20261110,
            'spearman_difference':rho(a,y)-rho(b,y),
            'descriptive_percentile_95':np.quantile(delta,[.025,.975]).tolist(),
            'caveat':'One protein and one mutant pair, correlated residues and gaps; interval is not independent biological replication.'}


def run():
    started=time.monotonic()
    assert sha(ROOT/'protocol.json')==PROTOCOL_SHA,'Frozen protocol differs'
    original=json.loads((ROOT.parent/'results/summary.json').read_text())
    for name,d in original['data'].items():assert sha(ROOT.parent/'data'/name)==d['sha256']
    scores=source_scores(ROOT.parent/'data/cypa-6BTA.cif')
    source_score_bytes=json.dumps(scores,sort_keys=True).encode()
    seq,a,ma,_=read_star(ROOT.parent/'data/cypa-27217.str')
    seq2,b,mb,_=read_star(ROOT.parent/'data/cypa-27218.str')
    assert [(i+1,x,z) for i,(x,z) in enumerate(zip(seq,seq2)) if x!=z]==[(97,'I','V')]
    all_rows=[];all_metrics={};uncertainty={};top=[]
    for temp in (283.,298.):
        rows=[]
        for pos,values in scores.items():
            if any((pos,atom) not in a[temp] or (pos,atom) not in b[temp] for atom in ('H','N')):continue
            dh=b[temp][pos,'H']-a[temp][pos,'H'];dn=b[temp][pos,'N']-a[temp][pos,'N']
            rows.append({'temperature_K':int(temp),'residue':pos,'csp':math.hypot(dh,dn/5),
                         'abs_delta_H':abs(dh),'abs_delta_N':abs(dn),**values})
        rows.sort(key=lambda r:r['residue'])
        all_rows+=rows
        for radius in (0,2,5):
            eligible=[r for r in rows if radius==0 or abs(r['residue']-97)>radius]
            all_metrics[f'{int(temp)}K_exclude_{radius}']={target:{name:metric(eligible,name,target) for name in SCORES}
                                                          for target in ('csp','abs_delta_H','abs_delta_N')}
            if radius==2:
                uncertainty[str(int(temp))]=bootstrap(eligible)
                measured=sorted(eligible,key=lambda r:-r['csp'])[:10]
                top.append({'temperature_K':int(temp),'largest_measured_csp_residues':[r['residue'] for r in measured]})
        # Cross-check CA baseline against original before evaluating improvement.
        old=original['metrics'][f'{int(temp)}K_exclude_2']['csp']['geometry_locality']
        new=all_metrics[f'{int(temp)}K_exclude_2']['csp']['geometry_locality']
        assert abs(old['spearman']-new['spearman'])<1e-12 and old['n']==new['n']
    carbon=[]
    for pos,values in scores.items():
        if abs(pos-97)<=2 or any((pos,atom) not in a[298.] or (pos,atom) not in b[298.] for atom in ('CA','CB')):continue
        dc=[b[298.][pos,atom]-a[298.][pos,atom] for atom in ('CA','CB')]
        carbon.append({'temperature_K':298,'residue':pos,'csp':math.hypot(*dc),**values})
    carbon.sort(key=lambda r:r['residue'])
    summary={'protocol_sha256':PROTOCOL_SHA,'source_structure_sha256':sha(ROOT.parent/'data/cypa-6BTA.cif'),
             'source_scores_sha256':hashlib.sha256(source_score_bytes).hexdigest(),
             'primary':PRIMARY,'metrics':all_metrics,'paired_block_bootstrap':uncertainty,
             'carbon_298K_secondary':{name:metric(carbon,name) for name in SCORES},
             'top_measured_sites':top,'elapsed_seconds':time.monotonic()-started,
             'point_improvement_both_temperatures':all(x['spearman_difference']>0 for x in uncertainty.values())}
    for name,rows in [('scores.csv',all_rows),('carbon_scores.csv',carbon)]:
        with (ROOT/name).open('w',newline='') as fh:
            w=csv.DictWriter(fh,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    (ROOT/'results.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'primary':{k:v['csp'] for k,v in all_metrics.items() if k.endswith('exclude_2')},
                      'bootstrap':uncertainty,'carbon':summary['carbon_298K_secondary'],
                      'point_improvement_both_temperatures':summary['point_improvement_both_temperatures']},indent=2))
    return summary


if __name__=='__main__':run()
