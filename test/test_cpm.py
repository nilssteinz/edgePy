import os
import unittest
import pandas as pd
from pandas import testing

from edgePy.cpm import CPM

DIR_PATH = os.path.split(os.getcwd())[0]
if "edgePy" not in DIR_PATH:
    DIR_PATH += "/edgePy"


class test_norm(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.main_data = pd.read_csv(DIR_PATH + "/data/data.txt", index_col=0, sep="\t")
        cls.main_factor = pd.read_csv(
            DIR_PATH + "/data/factors.txt", dtype=float, header=None
        )

    def setUp(self) -> None:
        self.data = self.main_data.copy()
        self.factor = self.main_factor.copy()
        self.Norm = CPM()

    def tearDown(self) -> None:
        pass

    def test_set_data(self):
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
        self.Norm.set_factor(self.factor)
        testing.assert_frame_equal(self.Norm.factor, self.factor)

    def test_set_factor_fail(self):
        with self.assertRaises(Exception) as e:
            self.Norm.set_factor(0)
        self.assertIsInstance(e.exception, (TypeError,))

    def test_CPM_not_norm(self):
        data = self.data.copy()
        result = CPM().cpm(data, normalized_lib_size=False)
        target = 2.751621e01
        self.assertEqual(target, result.iloc[0, 0], "CPM(CPM=F) calcualtions are off")

    def test_CPM_not_norm_log(self):
        data = self.data.copy()
        result = self.Norm.cpm(data, normalized_lib_size=False, log=True)
        target = 3.08325
        self.assertEqual(
            target, result.iloc[1, 1], "CPM(CPM=F, log=T) calcualtions are off"
        )

    def test_CPM_norm(self):
        data = self.data.copy()
        self.Norm.set_factor(self.factor)
        result = self.Norm.cpm(data, normalized_lib_size=True)
        target = 2.644534e01
        self.assertEqual(target, result.iloc[0, 0], "CPM calcualtions are off")

    def test_rpkm_non_norm(self):
        data = self.data.copy()
        result = self.Norm.rpkm(data, gene_length=100)
        target = 2.7516207e02
        self.assertEqual(target, result.iloc[0, 0], "RMPK calcualtions are off")

    def test_rpkm_not_norm_log(self):
        data = self.data.copy()
        self.Norm.round_results = 7
        result = self.Norm.rpkm(
            data, normalized_lib_size=False, log=True, gene_length=100
        )
        target = 8.1144958
        self.assertEqual(
            target, result.iloc[0, 0], "rpkm(CPM=F, log=T) calcualtions are off"
        )
        self.Norm.round_results = 5

    def test_rpkm_norm(self):
        data = self.data.copy()
        self.Norm.set_factor(self.factor)
        result = self.Norm.rpkm(data, normalized_lib_size=True, gene_length=100)
        target = 2.6445338e02
        self.assertEqual(target, result.iloc[0, 0], "rpkm(CPM=T) calcualtions are off")

    def test_rpkm_fail(self):
        with self.assertRaises(Exception) as e:
            self.Norm.rpkm(gene_length=None)
        self.assertIsInstance(e.exception, (TypeError,))

    def test_usage_self_data_cpm(self):
        self.Norm__ = CPM(self.data)
        self.Norm__.cpm()
        target = 2.751621e01
        self.assertEqual(target, self.Norm__.cpm_result.iloc[0, 0])

    def test_usage_self_data_rpkm(self):
        self.Norm__ = CPM(self.data)
        self.Norm__.rpkm(gene_length=100)
        target = 2.7516207e02
        self.assertEqual(target, self.Norm__.rpkm_result.iloc[0, 0])

    def test_rpkm_norm_log(self):
        data = self.data.copy()
        self.Norm.set_factor(self.factor)
        result = self.Norm.rpkm(
            data, normalized_lib_size=True, gene_length=100, log=True
        )
        target = 8.05757
        self.assertEqual(
            target, result.iloc[0, 0], "rpkm(CPM=T, log=T) calcualtions are off"
        )

    def test_cpm_raises_factor_error(self):
        with self.assertRaises(Exception) as e:
            self.Norm.cpm(normalized_lib_size=True)
        self.assertIsInstance(e.exception, (NotImplementedError,))

    def test_rpkm_raises_factor_error(self):
        with self.assertRaises(Exception) as e:
            self.Norm.rpkm(normalized_lib_size=True)
        self.assertIsInstance(e.exception, (NotImplementedError,))


if __name__ == "__main__":
    unittest.main()
