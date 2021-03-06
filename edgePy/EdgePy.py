"""
The EdgePy controller.

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

from __future__ import annotations
import copy
import pandas as pd

from cpm import CPM
from NormFactor import NormFactor


class EdgePy(object):
    """EdgePy main controller.

    -----------------------



    """

    CPM = CPM()
    NormFactor = NormFactor()

    def __init__(self, data: pd.DataFrame, classes: pd.Index = None):
        """
        Beginning of an EdgePy analysis.

        Parameters
        ----------
        data : pd.DataFrame
            This is the raw data from an RNA-seq experiment.
        classes: pd.Index
            The

        """
        self.data = copy.deepcopy(data)
        self.norm_data = None
        self.factor = pd.DataFrame([1] * len(data.columns))
        print(self.factor)

    def do_cpm(
        self,
        normalized_lib_size: bool = False,
        lib_size: pd.Index or slice = None,
        log: bool = False,
        prior_count: float = 2,
    ) -> EdgePy:
        """Execute the CPM for the data in the controller."""
        self.CPM.set_data(self.data)
        self.CPM.set_factor(self.factor)
        self.norm_data = self.CPM.cpm(
            normalized_lib_size=normalized_lib_size,
            lib_size=lib_size,
            log=log,
            prior_count=prior_count,
        )
        return self

    def do_rpkm(
        self,
        normalized_lib_size: bool = False,
        lib_size: pd.Index or slice = None,
        log: bool = False,
        prior_count: float = 2,
        gene_length: pd.DataFrame = None,
    ) -> EdgePy:
        """Execute the RPKM for the data in the controller."""
        self.CPM.set_data(self.data)
        self.CPM.set_factor(self.factor)
        self.norm_data = self.CPM.rpkm(
            normalized_lib_size=normalized_lib_size,
            lib_size=lib_size,
            log=log,
            prior_count=prior_count,
            gene_length=gene_length,
        )
        return self


if __name__ == "__main__":
    data = pd.DataFrame([[100, 100], [50, 150], [100, 50]])
    gene_lenght = [
        1000,
        500,
        1000,
        1000,
        500,
        1000,
        1000,
        500,
        1000,
        1000,
        500,
        1000,
        1000,
        500,
        1000,
    ]
    edge = EdgePy(data)
    print(edge.do_cpm())
    print(edge.do_rpkm(gene_length=pd.DataFrame([100])).norm_data)
    print(edge.CPM.cpm(data))
