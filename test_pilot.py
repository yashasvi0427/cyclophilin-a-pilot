import unittest
from pathlib import Path

import numpy as np

from pilot import ROOT, metric, read_star, read_structure


class MeasurementChecks(unittest.TestCase):
    def test_deposit_identity_and_conditions(self):
        double, a, ma, _ = read_star(ROOT / 'data/cypa-27217.str')
        triple, b, mb, _ = read_star(ROOT / 'data/cypa-27218.str')
        self.assertEqual([(i+1, x, y) for i, (x,y) in enumerate(zip(double,triple)) if x != y], [(97,'I','V')])
        self.assertEqual(set(a), {283.,298.})
        # Independently checked against the deposited atom rows, not generated outputs.
        self.assertEqual(a[283.][2,'H'], 8.528)
        self.assertEqual(b[283.][2,'H'], 8.507)
        self.assertEqual(ma[283.]['conditions']['pH'], 6.5)
        self.assertEqual(ma[283.]['conditions'],mb[283.]['conditions'])
        self.assertNotEqual(a[298.][3,'H'], a[283.][3,'H'])

    def test_target_structure_cannot_pass_as_source(self):
        double,*_ = read_star(ROOT / 'data/cypa-27217.str')
        with self.assertRaises(AssertionError):
            read_structure(ROOT / 'data/cypa-5WC7.cif', double)

    def test_source_mapping_and_mutation_atom(self):
        double,*_ = read_star(ROOT / 'data/cypa-27217.str')
        xyz, audit = read_structure(ROOT / 'data/cypa-6BTA.cif', double)
        self.assertEqual(audit['network_sites']['97']['aa'],'I')
        self.assertEqual(float(np.linalg.norm(xyz[97]-xyz[97])),0.)
        self.assertEqual(audit['temperature_K'],['273'])

    def test_known_rank_metrics(self):
        rows=[{'residue':i,'score':float(i),'csp':float(i)} for i in range(8)]
        m=metric(rows,'score')
        self.assertAlmostEqual(m['spearman'],1.)
        self.assertEqual(m['top_quartile_precision'],1.)
        self.assertEqual(m['chance_precision'],.25)
        for r in rows:r['score']=7-r['residue']
        m=metric(rows,'score')
        self.assertAlmostEqual(m['spearman'],-1.)
        self.assertEqual(m['top_quartile_precision'],0.)

    def test_constant_scores_and_tie_determinism(self):
        rows=[{'residue':i,'score':0.,'csp':float(i)} for i in range(8)]
        self.assertIsNone(metric(rows,'score')['spearman'])
        self.assertEqual(metric(rows,'score'),metric(list(reversed(rows)),'score'))


if __name__=='__main__':unittest.main()
