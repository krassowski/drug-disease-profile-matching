import numpy as np
from numba import jit


def identity(x):
    return x


@jit
def split_to_pos_and_neg(values):
    """Return indices of positive and negative values"""
    i = 0
    u = 0
    d = 0
    up_i = np.empty(values.size, dtype=np.int32)
    down_i = np.empty(values.size, dtype=np.int32)
    for vi in range(values.shape[0]):
        v = values[vi]
        if v > 0:
            up_i[u] = i
            u += 1
        elif v != 0:
            down_i[d] = i
            d += 1
        i += 1

    return up_i[:u], down_i[:d]


def one(x: list):
    assert len(x) == 1
    return x[0]
