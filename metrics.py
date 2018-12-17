from numpy import mean, sign, nan
from statistics import stdev
# for population/biased standard deviation use:
# from numpy import std as stdev

from numba import jit


@jit
def signal_to_noise(case, control):
    """Calculates SNR as ratio of means difference and deviation sum.

    Case and control has to be tuples or other hashable iterable.

    Assumes that there are:
        - at least two samples in both case and control
        - the samples have non-zero variation
    """
    dev_sum = (stdev(case) + stdev(control))
    return (
        (mean(case) - mean(control))
        /
        dev_sum
     ) if dev_sum else nan


@jit()
def fold_change(case, control):
    return mean(case) / mean(control)


@jit()
def signed_fold_change(case, control):
    return abs(mean(case) / mean(control)) * sign(mean(case))
