import unittest
import pandas as pd
from pandas import testing
from edgePy.NormFactor import NormFactor


class test_edgePy(unittest.TestCase):
    def setUp(self) -> None:

        self.data = pd.DataFrame([[1, 1]] * 10, dtype=float)
        self.data.iloc[-1, -1] = 100
        self.nromFact = NormFactor()
        pass

    def tearDown(self) -> None:
        pass

    def test_TMM(self):
        target = pd.DataFrame([[30289.13, 30289.13]] * 10)
        target.iloc[-1, -1] = 3028912.66

        # testing.assert_frame_equal(target, self.nromFact.calcNormFactor(self.data, method="TMM", ))


if __name__ == "__main__":
    unittest.main()
