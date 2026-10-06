import unittest
import numpy as np
from grid_analysis import allow_fine
from sequence_analysis import interpolation_map
from strict_analysis import matches_input,InvalidComparison
from run_grid import plan
class GridChecks(unittest.TestCase):
    def test_headless_capacity(self):
        from experiment import Config, DesktopConfig
        with self.assertRaises(ValueError):
            DesktopConfig(length_mm=150.,cells_across=160).validate()
        self.assertEqual(Config(length_mm=150.,cells_across=160).validate().cells_across,160)
        with self.assertRaises(ValueError):
            Config(length_mm=150.,cells_across=256).validate()
        with self.assertRaises(ValueError):
            Config(length_mm=150.,cells_across=160,time_scale=0.).validate()
    def test_gate(self):
        self.assertFalse(allow_fine([]))
        self.assertFalse(allow_fine([{'pair':'severe_N080__severe_N120','passed':False}]))
        self.assertFalse(allow_fine([{'pair':'severe_N080__severe_N120','passed':True,'smoke_only':True}]))
        self.assertTrue(allow_fine([{'pair':'severe_N080__severe_N120','passed':True}]))
    def test_affine_interpolation_and_solid_exclusion(self):
        a=dict(x_m=np.array([.2,.6]),y_m=np.array([.2,.6]),solid=np.zeros((2,2),bool))
        b=dict(x_m=np.array([0.,.5,1.]),y_m=np.array([0.,.5,1.]),solid=np.zeros((3,3),bool))
        x,y=np.meshgrid(b['x_m'],b['y_m']);X,mask,interp=interpolation_map(a,b)
        xx,yy=np.meshgrid(a['x_m'],a['y_m'])
        np.testing.assert_allclose(interp(2*x+3*y),2*xx+3*yy)
        b['solid'][0,0]=True
        _,mask,_=interpolation_map(a,b);self.assertFalse(mask[0,0]);self.assertTrue(mask[1,1])
    def test_wave_speed(self):
        cases=[s for s in plan() if s['kind']=='severe']
        speeds=[.004/s['cells']/s['dt_s'] for s in cases]
        np.testing.assert_allclose(speeds,[speeds[0]]*3,rtol=1e-12)
if __name__=='__main__':unittest.main()
