import copy

import numpy as np
import pandas as pd


class norm:
    """

    """

    round_value = 5

    def __init__(
        self, data: pd.DataFrame = None, factor: list or pd.DataFrame = None
    ) -> None:
        self.data = data
        self.factor = factor
        self.run: bool = False
        self.rpkm_result = None
        self.cpm_result = None

    def set_data(self, data: pd.DataFrame) -> None:
        """

        :param data:
        :return:
        """
        if type(data) != pd.DataFrame:
            raise TypeError("not a DataFrame used")
        self.data = data

    def set_factor(self, factor: list or pd.DataFrame = None) -> None:
        self.factor = factor

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
        lib_size: pd.DataFrame = data.sum(axis=0)
        if normalized_lib_size:
            lib_size = self.factor.to_numpy()[0] * lib_size
        adjusted_prior_count = (
            (int(log) * prior_count) * len(data.columns) * lib_size / lib_size.sum()
        )
        adjusted_lib_size = lib_size + 2 * adjusted_prior_count
        data: pd.DataFrame = (data + int(log) * adjusted_prior_count) / (
            adjusted_lib_size
        ) * 1000000
        if log:
            return np.log2(data).__round__(self.round_value)
        return data.__round__(self.round_value)

    def __rpkm(
        self,
        data: pd.DataFrame = None,
        log: bool = False,
        gene_length: list or pd.DataFrame = None,
        prior_count: int = 2,
    ) -> pd.DataFrame:
        """

        :param prior_count:
        :param data:
        :param log:
        :param gene_length:
        :param prior_count:
        :return:
        """
        k_b = gene_length / 1000
        data = data / (k_b)
        if log:
            return np.log2(data)
        return data

    def __lib_select(
        self, data: pd.DataFrame = None, lib_size: pd.Index or slice = None
    ) -> pd.DataFrame:
        """

        :param data:
        :param lib_size:
        :return:
        """
        data_cp = copy.deepcopy(data)
        if lib_size:
            data_cp = data_cp.iloc[lib_size]
        if len(data_cp.index) < 2:
            raise ArithmeticError("data set is to small for Normalisation")
        return data_cp

    def cpm(
        self,
        data: pd.DataFrame = None,
        normalized_lib_size: bool = False,
        lib_size: pd.Index or slice = None,
        log: bool = False,
        prior_count: float = 2,
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
        data_cp = self.__lib_select(data, lib_size)
        self.cpm_result = self.__cpm(data_cp, normalized_lib_size, log, prior_count)
        self.run = True
        return self.cpm_result

    def rpkm(
        self,
        data: pd.DataFrame = None,
        normalized_lib_size: bool = False,
        lib_size: pd.Index or slice = None,
        log: bool = False,
        prior_count: float = 2,
        gene_length: int = None,
    ) -> pd.DataFrame:
        """

        :param data:
        :param normalized_lib_size:
        :param lib_size:
        :param log:
        :param prior_count:
        :param gene_length:
        :return:
        """
        if gene_length == None:
            raise UserWarning("gene_length empty. must have a list or DataFrame")
        data_cp = self.__lib_select(data, lib_size)
        data_cp = self.__cpm(data_cp, normalized_lib_size, log, prior_count=prior_count)
        self.rpkm_result = self.__rpkm(data_cp, log=log, gene_length=gene_length)
        self.run = True
        return self.rpkm_result
