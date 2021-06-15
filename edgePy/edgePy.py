import copy

import pandas as pd

from cpm import cpm


class edgePy(object):
    norm = cpm()

    def __init__(self, data: pd.DataFrame):
        self.data = copy.deepcopy(data)
        self.norm_data = None
        self.factor = [1] * len(data.columns)
        print(self.factor)

    def do_cpm(
        self,
        normalized_lib_size: bool = False,
        lib_size: pd.Index or slice = None,
        log: bool = False,
        prior_count: float = 2,
    ) -> None:
        self.norm.set_data(self.data)
        self.norm.set_factor(self.factor)
        self.norm_data = self.norm.cpm(
            normalized_lib_size=normalized_lib_size,
            lib_size=lib_size,
            log=log,
            prior_count=prior_count,
        )

    def do_rpkm(
        self,
        normalized_lib_size: bool = False,
        lib_size: pd.Index or slice = None,
        log: bool = False,
        prior_count: float = 2,
        gene_length: int = None,
    ) -> None:
        self.norm.set_data(self.data)
        self.norm.set_factor(self.factor)
        self.norm_data = self.norm.rpkm(
            normalized_lib_size=normalized_lib_size,
            lib_size=lib_size,
            log=log,
            prior_count=prior_count,
            gene_length=gene_length,
        )

    def __str__(self):
        return data.__str__()


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
    edge = edgePy(data)
    print(edge.do_cpm())
    print(edge.do_rpkm(gene_length=100))
