from rpy2.rinterface._rinterface import RRuntimeError
from rpy2.robjects.packages import importr, data
from rpy2.robjects import Formula, r
from rpy2.robjects import StrVector, Vector
import re
from rpy2.robjects import pandas2ri
from rpy2.robjects.pandas2ri import ri2py as r2p
from rpy2.robjects.pandas2ri import py2ri as p2r


pandas2ri.activate()

source = r.source

base = importr('base')
c = base.c


def formula(f_str, df=None, **variables):
    assert not (df is not None and variables)
    fm = Formula(f_str)
    env = fm.environment
    if df is not None:
        variables = {
            name: df[name]
            for name in re.split('\W+', f_str)
        }
    for name, var in variables.items():
        if name not in f_str:
            raise ValueError(f'Variable {name} not used in {f_str}')
        env[name] = var
    return fm


class G:
    def __init__(self, plot):
        ggplot2 = importr("ggplot2")
        self.plot = plot
        self.add = ggplot2._env['%+%']

    def __add__(self, other):
        return G(self.add(self.plot, other))


def data_frame_from_matrix(matrix):
    from pandas import DataFrame
    import numpy as np
    
    rows_names = list(matrix.names[0])
    columns_names = list(matrix.names[1])

    return DataFrame(
        data=np.array(matrix),
        columns=columns_names,
        index=rows_names
    )


def one_or_all(x):
    if len(x) == 1:
        return x[0]
    return x


def r_ks_test(x, y, **kwargs):
    if len(x) < 2 or len(y) < 2:
        return {'p.value': 1}
    from pandas import Series
    try:
        result = r['ks.test'](p2r(Series(x)), p2r(Series(y)), **kwargs)
    except RRuntimeError:
        print(x)
        print(y)
        print(kwargs)
        raise
    return {
        key: one_or_all(value)
        for key, value in result.items()
    }
