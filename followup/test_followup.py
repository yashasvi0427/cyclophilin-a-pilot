import unittest
from pathlib import Path
import numpy as np
from run import minimum_distance,diffuse,source_scores,ROOT


class FollowupChecks(unittest.TestCase):
    def test_atom_distance_translation_rotation(self):
        a=np.array([[0.,0,0],[0,2,0]]);b=np.array([[3.,0,0]])
        self.assertAlmostEqual(minimum_distance(a,b),3)
        self.assertAlmostEqual(minimum_distance(a+19,b+19),3)
        q=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        self.assertAlmostEqual(minimum_distance(a@q,b@q),3)

    def test_diffusion_permutation_and_isolated_node(self):
        ca=np.array([[0.,0,0],[4,0,0],[8,0,0],[30,0,0]])
        p=diffuse(ca,0)
        self.assertAlmostEqual(p.sum(),1)
        self.assertEqual(p[3],0)
        self.assertTrue(np.all(p>=0))
        order=np.array([2,0,3,1])
        reordered=diffuse(ca[order],1)
        np.testing.assert_allclose(reordered,p[order],atol=1e-12)

    def test_rejects_post_mutation_geometry(self):
        with self.assertRaises(AssertionError):
            source_scores(ROOT.parent/'data/cypa-5WC7.cif')

    def test_source_predictor_has_no_label_inputs(self):
        values=source_scores(ROOT.parent/'data/cypa-6BTA.cif')
        self.assertEqual(values[97]['geometry_locality'],0)
        self.assertEqual(set(values[113]),{'geometry_locality','sidechain_amide_distance',
                                         'sidechain_ca_distance','all_heavy_distance','contact_diffusion'})
        self.assertTrue(all(np.isfinite(v) for r in values.values() for v in r.values()))


if __name__=='__main__':unittest.main()
