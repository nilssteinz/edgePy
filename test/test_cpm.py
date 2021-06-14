import os
import unittest

import numpy as np
import pandas as pd
from pandas import testing

from edgePy.cpm import norm

DIR_PATH = os.path.split(os.getcwd())[0]

class test_norm(unittest.TestCase):
    def setUp(self) -> None:

        self.main_data = pd.read_csv(DIR_PATH+"/data/data.txt", index_col=0, sep="\t")
        self.main_factor = pd.read_csv(DIR_PATH+"/data/factors.txt", dtype=float, header=None)
        self.Norm = norm()
        self.Norm.set_factor(self.main_factor)

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

    def test_CPM_not_norm(self):
        data = self.main_data.copy()
        result = norm().cpm(data, normalized_lib_size=False)
        target = 2.751621e01
        self.assertEqual(target, result.iloc[0, 0], "CPM(norm=F) calcualtions are off")

    def test_CPM_not_norm_log(self):
        data = self.main_data.copy()
        result = self.Norm.cpm(data, normalized_lib_size=False, log=True)
        target = 3.08325
        self.assertEqual(
            target, result.iloc[1, 1], "CPM(norm=F, log=T) calcualtions are off"
        )

    def test_CPM_norm(self):
        data = self.main_data.copy()
        result = self.Norm.cpm(data, normalized_lib_size=True)
        target = 2.644534e01
        self.assertEqual(target, result.iloc[0, 0], "CPM calcualtions are off")

    def test_rpkm_non_norm(self):
        data = self.main_data.copy()
        result = self.Norm.rpkm(data, gene_length=100)
        target = 2.7516207e02
        self.assertEqual(target, result.iloc[0, 0], "RMPK calcualtions are off")

    def test_rpkm_not_norm_log(self):
        data = self.main_data.copy()
        self.Norm.round_value = 7
        result = self.Norm.rpkm(data, normalized_lib_size=False, log=True, gene_length=100)
        target = 8.1144958
        self.assertEqual(
            target, result.iloc[0, 0], "rpkm(norm=F, log=T) calcualtions are off"
        )
        self.Norm.round_value = 5

    def test_rpkm_norm(self):
        data = self.main_data.copy()
        result = self.Norm.rpkm(
            data, normalized_lib_size=True, gene_length=100
        )
        target = 2.6445338e02
        self.assertEqual(target, result.iloc[0, 0], "rpkm(norm=T) calcualtions are off")


if __name__ == "__main__":
    unittest.main()
