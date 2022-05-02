#TODO
#test coordinate when calculate coordinate
#tested library score
#tested assign method (n2g, g2n)
import unittest
from . import grid as grid

class TestGenieGrid(unittest.TestCase):
    def test_n2g(self):
        self.assertEqual(grid.lon_n2g(175), -185)

    def test_g2n(self):
        self.assertEqual(grid.lon_g2n(-260), 100)
        self.assertEqual(grid.lon_g2n(-360), )

if __name__ == '__main__':
    unittest.main()
#method1: cal_rmse(xr)
#method2: cal_rmse(nc4)
#method3: xe.rmse