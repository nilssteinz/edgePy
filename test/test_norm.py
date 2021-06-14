import unittest

import numpy as np
import pandas as pd
from pandas import testing

from edgePy.norm import norm


class test_norm(unittest.TestCase):
    def setUp(self) -> None:
        self.main_data = pd.read_csv("../data/data.txt", index_col=0, sep="\t")
        self.Norm = norm()

    def tearDown(self) -> None:
        pass

    def test_set_data(self):
        """


        """
        self.assertEqual(self.Norm.data, None)
        data = pd.DataFrame([[10, 10], [10, 10]])
        self.Norm.set_data(data)
        testing.assert_frame_equal(self.Norm.data, data)

    def test_lib_size_cut_off_correct(self):
        data = pd.DataFrame(np.arange(0, 10))
        target_4_6 = norm().cpm(data.iloc[4:6])
        testing.assert_frame_equal(
            norm().cpm(data, lib_size=slice(4, 6)), target_4_6, check_dtype=False
        )

    def test_lib_size_raise_error_if_empty(self):
        """
        test if the __lib_select functions raises a error if you have less than 2 indexes.
        CPM should not return a thing if you only have a single row.
        :return:
        """
        data = pd.DataFrame(np.arange(0, 10))
        with self.assertRaises(Exception) as e:
            norm().cpm(data, lib_size=slice(7, 2))
        self.assertTrue(type(e.exception) in [ArithmeticError])


    def test_CPM_non_norm(self):
        data = self.main_data.copy()
        result = norm().cpm(data)
        target = 2.751621e+01
        self.assertEqual(target,result.iloc[0,0], "CPM calcualtions are off")
        pass

    def test_rpkm_non_norm(self):
        data = self.main_data.copy()
        result = norm().rpkm(data, gene_length=100)
        target = 2.751621e+02
        self.assertEqual(target,result.iloc[0,0], "RMPK calcualtions are off")
        pass


if __name__ == "__main__":
    unittest.main()
