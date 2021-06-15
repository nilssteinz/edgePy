"""python version of the edgeR package CPM <Counts Per Million>.

Copyright (C) 2021 Nils Steinz <nils.steinz@hotmail.com>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 License as published by the Free Software Foundation, either version 3 of the License,
  or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.
 If not, see <https://www.gnu.org/licenses/>.
"""


import copy

import numpy as np
import pandas as pd


class CPM:
    """
    python version of the edgeR package CPM <Counts Per Million>.
    -------------------------------------------------------------


    

    parameters
    ----------
    round_results: int
    """

    round_results: int = 5

    def __init__(
        self, data: pd.DataFrame = None, factor: list or pd.DataFrame = None
    ) -> None:
        """

        :param data:
        :param factor:
        """
        self.data = data
        self.factor = factor
        self.rpkm_result = None
        self.cpm_result = None

    def set_data(self, data: pd.DataFrame) -> None:
        """

        :param data:
        :return:
        """
        if not isinstance(data, pd.DataFrame):
            raise TypeError("not a DataFrame used")
        self.data = data

    def set_factor(self, factor: list or pd.DataFrame) -> None:
        """

        :param factor:
        :return:
        """
        if not isinstance(factor, pd.DataFrame):
            raise TypeError("not a DataFrame used")
        self.factor = factor

    def __calculations(
        self,
        data: pd.DataFrame = None,
        lib_size: pd.DataFrame = None,
        normalized_lib_size: bool = False,
        log: bool = False,
        prior_count: float = 2,
        gene_length: pd.Series or pd.DataFrame = 1000,
    ) -> pd.DataFrame:
        if lib_size is None:
            lib_size: pd.DataFrame = data.sum(axis=0)
        if normalized_lib_size:
            lib_size = self.factor.to_numpy()[0] * lib_size

        adjusted_prior_count = (
            (int(log) * prior_count) * len(data.columns) * lib_size / lib_size.sum()
        )
        adjusted_lib_size = lib_size + 2 * adjusted_prior_count
        data: pd.DataFrame = (
            (
                (data + int(log) * adjusted_prior_count) / (adjusted_lib_size) * 1000000
            ).transpose()
            / (gene_length / 1000)
        ).transpose()
        if log:
            return np.log2(data).__round__(self.round_results)
        return data.__round__(self.round_results)

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
        if normalized_lib_size and self.factor is None:
            raise NotImplementedError(
                "cannot use 'normalized_lib_size' without factor values."
            )
        if data is None:
            data = self.data
        data_cp = data.copy()
        self.cpm_result = self.__calculations(
            data=data_cp,
            lib_size=lib_size,
            normalized_lib_size=normalized_lib_size,
            log=log,
            prior_count=prior_count,
        )
        return self.cpm_result

    def rpkm(
        self,
        data: pd.DataFrame = None,
        normalized_lib_size: bool = False,
        lib_size: pd.Index or slice = None,
        log: bool = False,
        prior_count: float = 2,
        gene_length: list or pd.DataFrame = None,
    ) -> pd.DataFrame:
        """Reads per kilobase million



        Parameters
        ----------
        :param data: pd.Dataframe
            Read counts DataFrame.
        :param normalized_lib_size: boolean
            Use normalization factor to adjust the library size.
        :param lib_size: pd.DataFrame
            A Dataframe that contains the library size of the reads.
        :param log: boolean
        :param prior_count: float
            value added to remedy -inf with log. only functions when Log=True.
        :param gene_length: pd.DataFrame
            A list that contains the size of all the genes used in the data to
        :return: pd.DataFrame
            returns a DataFrame of same size as input where the reads per kilobase million has been calculated.

        """
        if normalized_lib_size and self.factor is None:
            raise NotImplementedError(
                "cannot use 'normalized_lib_size' without factor values."
            )
        elif isinstance(gene_length, (type(None),)):
            raise TypeError("gene_length empty. must have a list or DataFrame")
        elif data is None:
            data = self.data
        data_cp = data.copy()
        self.rpkm_result = self.__calculations(
            data_cp,
            log=log,
            gene_length=gene_length,
            prior_count=prior_count,
            normalized_lib_size=normalized_lib_size,
        )
        return self.rpkm_result
