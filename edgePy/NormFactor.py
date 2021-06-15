"""
Calculator for the normalisation factor used by edgeR <calcNormFactors.R>.

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
import pandas as pd
import numpy as np


class NormFactor:
    """Calculator for the normalisation factor used by edgeR <calcNormFactors.R>."""

    def __init__(self, data: pd.DataFrame = None) -> None:
        """Initialize setup for the method usage."""
        self.data = data

    def calcNormFactor(
        self,
        data: pd.DataFrame = None,
        method: str = "TMM",
        libSize: pd.DataFrame = 1,
        refColumn: pd.Index = None,
        logratioTrim: float = 0.3,
        sumTrim: float = 0.05,
        doWeighting: bool = True,
        Acutoff: float = -1e10,
        p: float = 0.75,
    ) -> pd.DataFrame:
        """Beginning of the Calculator."""
        data_dp = data.copy()

        return data_dp
