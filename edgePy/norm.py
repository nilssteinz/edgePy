import copy

import numpy as np
import pandas as pd


class norm:
    """

    """
    def __init__(self, data: pd.DataFrame = None):
        self.data = data

    def set_data(self, data: pd.DataFrame):
        """

        :param data:
        :return:
        """
        if type(data) != pd.DataFrame:
            raise TypeError("not a DataFrame used")
        self.data = data

    def _cpm(
        self,
        data: pd.DataFrame = None,
        normalized_lib_size: bool = False,
        log: bool = False,
        prior_count: float = 2,
    ):
        """

        :param data:
        :param normalized_lib_size:
        :param log:
        :param prior_count:
        :return:
        """
        lib_size = data.sum(axis=0)
        adjusted_lib_size = lib_size + 2 * (prior_count * lib_size / lib_size.mean())
        data = data / (adjusted_lib_size) * 100000
        if log:
            return np.log2(data)
        return data

    def _rmpk(
        self,
        data: pd.DataFrame = None,
        log: bool = False,
        gene_length: list or pd.DataFrame = None,
    ):
        """

        :param data:
        :param log:
        :param gene_length:
        :return:
        """
        for x in data.index:
            k_b = gene_length[x] / 1000
            data.loc[x, :] = data.loc[x, :] / (k_b)
        if log:
            return np.log(data)
        return data

    def __lib_size(self, data: pd.DataFrame = None, lib_size: pd.Index or slice = None):
        """

        :param data:
        :param lib_size:
        :return:
        """
        data_cp = copy.deepcopy(data)
        if lib_size:
            data_cp = data_cp.iloc[lib_size]
        return data_cp

    def cmp(
        self,
        data: pd.DataFrame = None,
        normalized_lib_size: bool = False,
        lib_size: pd.Index or slice = None,
        log: bool = False,
        prior_count: float = 0,
    ):
        """

        :param data:
        :param normalized_lib_size:
        :param lib_size:
        :param log:
        :param prior_count:
        :return:
        """
        if data is None:
            data = self.data
            print(data)
        data_cp = self.__lib_size(data, lib_size)
        data_cp = self._cpm(data_cp, normalized_lib_size, log, prior_count)

        return data_cp

    def rmpk(
        self,
        data: pd.DataFrame = None,
        normalized_lib_size: bool = False,
        lib_size: pd.Index or slice = None,
        log: bool = False,
        prior_count: float = 0,
        gene_length: list = None,
    ):
        data_cp = self.__lib_size(data, lib_size)
        data_cp = self._cpm(
            data_cp, normalized_lib_size, log=False, prior_count=prior_count
        )
        data = self._rmpk(data_cp, log, gene_length)

        return data

    def __call__(self, *args, **kwargs):
        if self.data:
            return self.cmp()
        else:
            raise Exception("wrong call")
