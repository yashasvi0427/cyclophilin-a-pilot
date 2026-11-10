"""Retrospective CypA measurement/localization pilot. No fitted predictor."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import time
import urllib.request
from pathlib import Path

import numpy as np
import pynmrstar
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
SOURCES = {
    **{f"cypa-{i}.str": f"https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr{i}/bmr{i}_3.str" for i in (27217, 27218)},
    **{f"cypa-{i}.cif": f"https://files.rcsb.org/download/{i}.cif" for i in ("3K0N", "3K0O", "6BTA", "5WC7")},
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fetch():
    (ROOT / "data").mkdir(exist_ok=True)
    for name, url in SOURCES.items():
        path = ROOT / "data" / name
        if not path.exists():
            urllib.request.urlretrieve(url, path)
    manifest = {name: {"url": url, "bytes": (ROOT / "data" / name).stat().st_size,
                       "sha256": sha(ROOT / "data" / name)} for name, url in SOURCES.items()}
    return manifest


def read_star(path):
    entry = pynmrstar.Entry.from_file(str(path))
    entity = entry.get_loops_by_category("Entity_comp_index")[0]
    seq_rows = entity.get_tag(["ID", "Comp_ID"])
    sequence = "".join(seq1(aa) for _, aa in sorted(seq_rows, key=lambda r: int(r[0])))
    conditions = {}
    for loop in entry.get_loops_by_category("Sample_condition_variable"):
        rows = loop.get_tag(["Type", "Val", "Sample_condition_list_ID"])
        for kind, val, ident in rows:
            conditions.setdefault(ident, {})[kind] = float(val)
    shifts, metadata, skipped = {}, {}, {}
    for frame in entry:
        if frame.category != "assigned_chemical_shifts":
            continue
        cond = frame.get_tag("_Assigned_chem_shift_list.Sample_condition_list_ID")[0]
        temp = conditions[cond]["temperature"]
        shifts[temp] = {}
        skipped[temp] = 0
        details = frame.get_tag("_Assigned_chem_shift_list.Details")[0]
        metadata[temp] = {"conditions": conditions[cond], "details": details}
        loop = frame.get_loop_by_category("Atom_chem_shift")
        for row in loop:
            data = dict(zip(loop.tags, row))
            if data["Ambiguity_code"] != "1" or data["Val"] in (".", "?"):
                skipped[temp] += 1
                continue
            pos, author = int(data["Seq_ID"]), int(data["Auth_seq_ID"])
            assert pos == author, "Author/NMR numbering mismatch"
            assert sequence[pos - 1] == seq1(data["Comp_ID"])
            key = pos, data["Atom_ID"]
            assert key not in shifts[temp], f"Duplicate assignment {key}"
            shifts[temp][key] = float(data["Val"])
    return sequence, shifts, metadata, skipped


def read_structure(path, expected_sequence=None):
    d = MMCIF2Dict(str(path))
    structure = MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    chain = next(structure[0].get_chains())
    coords, residues, alternates = {}, {}, {}
    for residue in chain:
        if residue.id[0] != " ":
            continue
        pos = residue.id[1]
        if "CA" not in residue:
            continue
        aa = seq1(residue.resname)
        if expected_sequence is not None:
            assert expected_sequence[pos - 1] == aa, f"PDB sequence mismatch {path}:{pos}"
        coords[pos] = residue["CA"].coord.astype(float)
        residues[pos] = aa
        alternatives = {}
        for atom in residue.get_unpacked_list():
            alt = atom.get_altloc().strip()
            if alt:
                alternatives.setdefault(alt, set()).add(atom.get_name())
        if alternatives:
            alternates[pos] = {key: sorted(v) for key, v in alternatives.items()}
    audit = {"chain": chain.id, "temperature_K": d.get("_diffrn.ambient_temp", []),
             "n_CA_residues": len(coords), "residues_with_alternate_atoms": len(alternates),
             "network_sites": {str(i): {"aa": residues.get(i), "alternate_atoms": alternates.get(i, {})}
                               for i in (97, 99, 113, 115)}}
    return coords, audit


def esm_scores(sequences, output, model_cache):
    import esm
    import torch
    torch.set_num_threads(2)
    torch.manual_seed(0)
    torch.hub.set_dir(str(model_cache))
    started = time.monotonic()
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    model.eval()
    # CPU chosen deliberately: short sequence, deterministic, modest unified memory use.
    _, _, tokens = alphabet.get_batch_converter()([("double", sequences[0]), ("triple", sequences[1])])
    with torch.inference_mode():
        embeddings = model(tokens, repr_layers=[6])["representations"][6][:, 1:-1].cpu().numpy()
    values = np.linalg.norm(embeddings[1] - embeddings[0], axis=1)
    np.savez_compressed(output / "esm_embeddings.npz", double=embeddings[0], triple=embeddings[1])
    ckpt = model_cache / "checkpoints" / "esm2_t6_8M_UR50D.pt"
    metadata = {"name": "esm2_t6_8M_UR50D", "layer": 6, "parameters": sum(p.numel() for p in model.parameters()),
                "device": "cpu", "threads": 2, "elapsed_seconds": time.monotonic() - started,
                "torch_version": torch.__version__, "weights_sha256": sha(ckpt),
                "weights_bytes": ckpt.stat().st_size,
                "source": "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t6_8M_UR50D.pt"}
    return values, metadata


def metric(rows, score, target="csp"):
    x, y = np.array([r[score] for r in rows]), np.array([r[target] for r in rows])
    n, k = len(rows), math.ceil(len(rows) / 4)
    # Explicit deterministic tie breaking instead of an arbitrary percentile threshold.
    truth = {r["residue"] for r in sorted(rows, key=lambda r: (-r[target], r["residue"]))[:k]}
    predicted = {r["residue"] for r in sorted(rows, key=lambda r: (-r[score], r["residue"]))[:k]}
    rho = None if np.ptp(x) == 0 or np.ptp(y) == 0 else float(spearmanr(x, y).statistic)
    return {"n": n, "k": k, "spearman": rho, "top_quartile_precision": len(truth & predicted) / k,
            "chance_precision": k / n, "score_ties": n - len(set(x.tolist()))}


def run(use_esm=False, cached_esm=False, model_cache=None):
    started = time.monotonic()
    output = ROOT / "results"
    output.mkdir(exist_ok=True)
    manifest = fetch()
    double, a, ma, sa = read_star(ROOT / "data/cypa-27217.str")
    triple, b, mb, sb = read_star(ROOT / "data/cypa-27218.str")
    changes = [(i + 1, x, y) for i, (x, y) in enumerate(zip(double, triple)) if x != y]
    assert len(double) == len(triple) == 165 and changes == [(97, "I", "V")], changes
    assert set(a) == set(b) == {283., 298.}
    for temp in a:
        assert ma[temp]["conditions"] == mb[temp]["conditions"]
        assert " ".join(ma[temp]["details"].split()) == " ".join(mb[temp]["details"].split())
    coords, source_audit = read_structure(ROOT / "data/cypa-6BTA.cif", double)
    _, target_audit = read_structure(ROOT / "data/cypa-5WC7.cif", triple)
    model_values = None
    model_meta = {"status": "not_run"}
    if use_esm:
        model_values, model_meta = esm_scores((double, triple), output, Path(model_cache))
    elif cached_esm:
        emb = np.load(output / "esm_embeddings.npz")
        model_values = np.linalg.norm(emb["triple"] - emb["double"], axis=1)
        model_meta = json.loads((output / "summary.json").read_text())["model"]
    scores = ["sequence_locality", "geometry_locality"] + (["esm2_embedding_change"] if model_values is not None else [])
    all_rows, results, matching = [], {}, {}
    for temp in sorted(a):
        matched, excluded = [], []
        for pos in range(1, 166):
            if pos not in coords or any((pos, atom) not in a[temp] or (pos, atom) not in b[temp] for atom in ("H", "N")):
                excluded.append(pos)
                continue
            dh = b[temp][pos, "H"] - a[temp][pos, "H"]
            dn = b[temp][pos, "N"] - a[temp][pos, "N"]
            row = {"temperature_K": temp, "residue": pos, "aa_double": double[pos-1], "aa_triple": triple[pos-1],
                   "delta_H_ppm": dh, "delta_N_ppm": dn, "abs_delta_H": abs(dh), "abs_delta_N": abs(dn),
                   "csp": math.sqrt(dh * dh + (dn / 5) ** 2),
                   "sequence_locality": -abs(pos - 97),
                   "geometry_locality": -float(np.linalg.norm(coords[pos] - coords[97]))}
            if model_values is not None:
                row["esm2_embedding_change"] = float(model_values[pos - 1])
            matched.append(row)
        all_rows.extend(matched)
        matching[str(int(temp))] = {"matched": len(matched), "excluded_missing": excluded,
                                   "skipped_ambiguous_assignments_double": sa[temp], "skipped_ambiguous_assignments_triple": sb[temp]}
        for exclusion in (0, 2, 5):
            # 0 denotes all matched sites, for sensitivity only.
            rows = [r for r in matched if exclusion == 0 or abs(r["residue"] - 97) > exclusion]
            results[f"{int(temp)}K_exclude_{exclusion}"] = {
                target: {score: metric(rows, score, target) for score in scores}
                for target in ("csp", "abs_delta_H", "abs_delta_N")}
    with (output / "residue_scores.csv").open("w") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    summary = {"protocol_sha256": sha(ROOT / "PROTOCOL.md"), "data": manifest,
               "sequences": {"double": double, "triple": triple, "changes": changes},
               "conditions": {"double": ma, "triple": mb}, "matching": matching,
               "structure_audit": {"source_6BTA": source_audit, "target_5WC7": target_audit},
               "model": model_meta, "metrics": results,
               "environment": {"python": platform.python_version(), "platform": platform.platform(),
                               "numpy": np.__version__}, "elapsed_seconds": time.monotonic() - started}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"matching": matching, "primary": {k: v["csp"] for k, v in results.items() if k.endswith("exclude_2")},
                      "model": model_meta}, indent=2))
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--esm", action="store_true")
    p.add_argument("--cached-esm", action="store_true")
    p.add_argument("--model-cache", default="/Volumes/Yashasvi/Research/cache/torch")
    args = p.parse_args()
    run(args.esm, args.cached_esm, args.model_cache)
