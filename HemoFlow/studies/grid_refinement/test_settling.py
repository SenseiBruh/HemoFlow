import unittest
import numpy as np
import strict_analysis as A
from settling_analysis import cycle_rows
from run_grid import plan

class Checks(unittest.TestCase):
    def data(self,drift=0.):
        t=np.linspace(0,4,4001);y=2+np.sin(2*np.pi*t)+drift*t
        d={'time_s':t,**{k:y.copy() for k in A.SCALARS}}
        p={'time_s':t,**{f'{w}_x{x}_pa':y.copy() for w in ('lower','upper') for x in (20,28,30,32,40,50,55)}}
        return d,p
    def test_plan(self):
        p=plan();self.assertEqual([s['key'] for s in p],['healthy_N080','healthy_N120','severe_N080','severe_N120','healthy_N160','severe_N160'])
        self.assertEqual(p[2]['cycles'],10)
        self.assertAlmostEqual(p[2]['cells']*p[2]['dt_s'],p[3]['cells']*p[3]['dt_s'])
    def test_periodic(self):
        _,rows=cycle_rows(*self.data(),1,4)
        self.assertEqual(len(rows),3)
        self.assertLess(max(v for r in rows for k,v in r.items() if k.endswith('_percent')),1e-10)
    def test_drift(self):
        _,rows=cycle_rows(*self.data(.2),1,4)
        self.assertGreater(rows[-1]['signed_wss_full_l2_percent'],1)
    def test_incomplete(self):
        with self.assertRaises(A.InvalidComparison):cycle_rows(*self.data(),1,5)
if __name__=='__main__':unittest.main()
