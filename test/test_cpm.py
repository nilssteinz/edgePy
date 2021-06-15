import os
import unittest

import numpy as np
import pandas as pd
from pandas import testing


from edgePy.cpm import cpm

DIR_PATH = os.path.split(os.getcwd())[0]
if "edgePy" not in DIR_PATH:
    DIR_PATH += "/edgePy"

class test_norm(unittest.TestCase):
    def setUp(self) -> None:

        self.main_data = pd.read_csv(DIR_PATH + "/data/data.txt", index_col=0, sep="\t")
        self.main_factor = pd.read_csv(
            DIR_PATH + "/data/factors.txt", dtype=float, header=None
        )
        self.Norm = cpm()

    def tearDown(self) -> None:
        pass

    def test_set_data(self):
        """

        """
        self.assertEqual(self.Norm.data, None)
        data = pd.DataFrame([[10, 10], [10, 10]])
        self.Norm.set_data(data)
        testing.assert_frame_equal(self.Norm.data, data)

    def test_set_data_fail(self):
        with self.assertRaises(Exception) as e:
            self.Norm.set_data(0)
        self.assertIsInstance(e.exception, (TypeError,))

    def test_set_factor(self):
        self.assertEqual(self.Norm.factor, None)
        data = pd.DataFrame([10, 10, 10, 10])
        self.Norm.set_factor(data)
        testing.assert_frame_equal(self.Norm.factor, data)

    def test_set_factor_fail(self):
        with self.assertRaises(Exception) as e:
            self.Norm.set_factor(0)
        self.assertIsInstance(e.exception, (TypeError,))

    def test_lib_size_cut_off_correct(self):
        data = pd.DataFrame(np.arange(0, 10))
        target_4_6 = self.Norm.cpm(data.iloc[4:6])
        testing.assert_frame_equal(
            cpm().cpm(data, lib_size=slice(4, 6)), target_4_6, check_dtype=False
        )

    def test_lib_size_raise_error_if_empty(self):
        data = pd.DataFrame(np.arange(0, 10))
        with self.assertRaises(Exception) as e:
            cpm().cpm(data, lib_size=slice(7, 2))
        self.assertIsInstance(e.exception, (ArithmeticError,))

    def test_CPM_not_norm(self):
        data = self.main_data.copy()
        result = cpm().cpm(data, normalized_lib_size=False)
        target = 2.751621e01
        self.assertEqual(target, result.iloc[0, 0], "CPM(cpm=F) calcualtions are off")

    def test_CPM_not_norm_log(self):
        data = self.main_data.copy()
        result = self.Norm.cpm(data, normalized_lib_size=False, log=True)
        target = 3.08325
        self.assertEqual(
            target, result.iloc[1, 1], "CPM(cpm=F, log=T) calcualtions are off"
        )

    def test_CPM_norm(self):
        data = self.main_data.copy()
        self.Norm.set_factor(self.main_factor)
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
        result = self.Norm.rpkm(
            data, normalized_lib_size=False, log=True, gene_length=100
        )
        target = 8.1144958
        self.assertEqual(
            target, result.iloc[0, 0], "rpkm(cpm=F, log=T) calcualtions are off"
        )
        self.Norm.round_value = 5

    def test_rpkm_norm(self):
        data = self.main_data.copy()
        self.Norm.set_factor(self.main_factor)
        result = self.Norm.rpkm(data, normalized_lib_size=True, gene_length=100)
        target = 2.6445338e02
        self.assertEqual(target, result.iloc[0, 0], "rpkm(cpm=T) calcualtions are off")

    def test_rpkm_fail(self):
        with self.assertRaises(Exception) as e:
            self.Norm.rpkm(gene_length=None)
        self.assertIsInstance(e.exception, (TypeError,))

    def test_usage_self_data_cpm(self):
        self.Norm__ = cpm(self.main_data)
        self.Norm__.cpm()
        result = self.Norm.cpm(self.main_data)
        testing.assert_frame_equal(self.Norm__.cpm_result, result)

    def test_usage_self_data_rpkm(self):
        self.Norm__ = cpm(self.main_data)
        self.Norm__.rpkm(gene_length=100)
        result = self.Norm.rpkm(data=self.main_data, gene_length=100)
        testing.assert_frame_equal(self.Norm__.rpkm_result, result)

    def test_rpkm_norm_log(self):
        data = self.main_data.copy()
        self.Norm.set_factor(self.main_factor)
        result = self.Norm.rpkm(
            data, normalized_lib_size=True, gene_length=100, log=True
        )
        target = 8.05757
        self.assertEqual(
            target, result.iloc[0, 0], "rpkm(cpm=T, log=T) calcualtions are off"
        )


if __name__ == "__main__":
    unittest.main()
