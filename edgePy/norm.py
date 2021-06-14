import copy

import numpy as np
import pandas as pd


class norm:
    """

    """
    def __init__(self, data: pd.DataFrame = None) -> None:
        self.data = data

    def set_data(self, data: pd.DataFrame) -> None:
        """

        :param data:
        :return:
        """
        if type(data) != pd.DataFrame:
            raise TypeError("not a DataFrame used")
        self.data = data

    def __cpm(
        self,
        data: pd.DataFrame = None,
        normalized_lib_size: bool = False,
        log: bool = False,
        prior_count: float = 2,
    ) -> pd.DataFrame:
        """

        :param data:
        :param normalized_lib_size:
        :param log:
        :param prior_count:
        :return: data
        """
        lib_size = data.sum(axis=0)
        adjusted_lib_size = lib_size + len(data.index) * (prior_count * lib_size / lib_size.sum())
        data = data / (adjusted_lib_size) * 1000000
        if log:
            return np.log2(data).__round__(5)
        return data.__round__(5)

    def __rpkm(
        self,
        data: pd.DataFrame = None,
        log: bool = False,
        gene_length: list or pd.DataFrame = None,
    ) -> pd.DataFrame:
        """

        :param data:
        :param log:
        :param gene_length:
        :return:
        """
        k_b = gene_length / 1000
        data = data / (k_b)
        if log:
            return np.log(data)
        return data

    def __lib_select(self, data: pd.DataFrame = None, lib_size: pd.Index or slice = None) -> pd.DataFrame:
        """

        :param data:
        :param lib_size:
        :return:
        """
        data_cp = copy.deepcopy(data)
        if lib_size:
            data_cp = data_cp.iloc[lib_size]
        if len(data_cp.index) <2:
            raise ArithmeticError("data set is to small for Normalisation")
        return data_cp

    def cpm(
        self,
        data: pd.DataFrame = None,
        normalized_lib_size: bool = False,
        lib_size: pd.Index or slice = None,
        log: bool = False,
        prior_count: float = 0,
    ) -> pd.DataFrame:
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
        data_cp = self.__lib_select(data, lib_size)
        data_cp = self.__cpm(data_cp, normalized_lib_size, log, prior_count)

        return data_cp

    def rpkm(
        self,
        data: pd.DataFrame = None,
        normalized_lib_size: bool = False,
        lib_size: pd.Index or slice = None,
        log: bool = False,
        prior_count: float = 0,
        gene_length: int = None,
    ) -> pd.DataFrame:
        if gene_length == None:
            raise UserWarning("gene_lenght empty. must have a list or DataFrame")
        data_cp = self.__lib_select(data, lib_size)
        data_cp = self.__cpm(
            data_cp, normalized_lib_size, log=False, prior_count=prior_count
        )
        data = self.__rpkm(data_cp, log, gene_length)

        return data

    def __call__(self, *args, **kwargs):
        if self.data:
            return self.cpm()
        else:
            raise Exception("wrong call")
