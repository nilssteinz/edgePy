import copy

import pandas as pd

from edgePy.norm import norm


class edgePy(object):

    norm = norm()

    def __init__(self, data: pd.DataFrame):
        self.data = copy.deepcopy(data)
        self.norm.set_data(data)

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
    print(edgePy(data).norm.cpm())
    print(edgePy.norm.cpm(data=data))
    print(edgePy.norm.rpkm(data=data, gene_length=[1000, 500, 1000]))
