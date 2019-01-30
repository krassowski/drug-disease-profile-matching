import numpy

from pandas import Series
from scipy.stats import combine_pvalues
from sklearn.metrics import roc_auc_score

from helpers import on_division_by_zero


def f1(precision, recall):
    return 2 * (precision * recall) / (precision + recall) if precision + recall else 0


@on_division_by_zero(fill_with=0)
def precision(true_positives, all_selected):
    return len(true_positives) / len(all_selected)


@on_division_by_zero(fill_with=0)
def recall(true_positives, all_positives):
    return len(true_positives) / len(all_positives)


def generalized_roc_auc_score(result):
    if len(set(result.expected)) < 2:
        return numpy.nan

    return roc_auc_score(result.expected, result.observed.tolist())


def fisher_method(pvalues):
    return combine_pvalues(pvalues, method='fisher')[1]


@on_division_by_zero(fill_with=numpy.inf)
def normalized_means_difference(a: Series, b: Series):
    _min = min([a.min(), b.min()])
    _max = min([a.max(), b.max()])
    scores_range = abs(_min - _max)
    return (a.mean() - b.mean()) / scores_range
