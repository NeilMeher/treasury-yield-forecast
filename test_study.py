import unittest
import numpy as np

from study import MATS, forecast, loadings, synthetic, walk


class ForecastTests(unittest.TestCase):
    def test_loading_endpoints_and_curvature(self):
        x=loadings([0,2.5,10000])
        np.testing.assert_allclose(x[0],[1,1,0])
        self.assertLess(x[-1,1],.001)
        self.assertGreater(x[1,2],x[-1,2])

    def test_unchanged_and_flat_curve(self):
        curves=np.full((80,len(MATS)),4.5)
        preds,factors=walk(curves,30)
        np.testing.assert_allclose(preds[12][0],4.5,atol=1e-11)
        np.testing.assert_allclose(factors[:,1:],0,atol=1e-11)

    def test_forecasts_unchanged_when_altering_future(self):
        _,a=synthetic()
        b=a.copy()
        b[61:]+=11
        left,_=walk(a,60)
        right,_=walk(b,60)
        np.testing.assert_allclose(left[12][0][:2],right[12][0][:2])
        self.assertFalse(np.allclose(left[12][1][:2],right[12][1][:2]))

    def test_ar_iteration(self):
        values=np.array([[8.],[5.],[3.5],[2.75]])
        self.assertAlmostEqual(forecast(values,2)[0],2.1875)
