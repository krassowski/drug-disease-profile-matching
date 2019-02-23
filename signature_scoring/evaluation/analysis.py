from pandas import DataFrame
from sklearn.metrics import log_loss

from helpers.gui import NeatNamespace
from helpers.r import r_ks_test

from .calculation_utilities import generalized_roc_auc_score


def ks_distance(a, b):
    return r_ks_test(a, b)['statistic']


def cluster_scores(scores, normalize=True, distance=ks_distance):
    functions = list(scores.func.unique())
    matrix = DataFrame(index=functions, columns=functions)
    for i, func in enumerate(functions):
        for other_func in functions[i+1:]:
            a = scores[scores.func == func]['score']
            b = scores[scores.func == other_func]['score']
            d = distance(a, b)
            matrix.loc[func, other_func] = d
            matrix.loc[other_func, func] = d
        matrix.loc[func, func] = 0
    if normalize:
        matrix = matrix / matrix.max().max()
    return matrix


def compute_metrics(known_status, expected='is_indication'):
    comparison = NeatNamespace(
        expected=list(known_status[expected]),
        observed=known_status.score
    )

    return {
        'auc': generalized_roc_auc_score(comparison),
        'log_loss': log_loss(known_status[expected].tolist(), known_status.score.tolist()),
        'auc0.01': generalized_roc_auc_score(comparison, max_fpr=0.01),
        'auc0.1': generalized_roc_auc_score(comparison, max_fpr=0.1)
    }


def compute_metrics_by_func(scores, without_unassigned):
    metrics_by_func = []
    for func in scores.func.unique():
        func_scores = scores[scores.func == func]
        if without_unassigned:
            known_scores = func_scores[func_scores.group != 'unassigned']
            func_scores = known_scores
        metrics_by_func.append(
            {'func': func, **compute_metrics(func_scores)}
        )
    return DataFrame(metrics_by_func)


def normalize_scores(scores, rescale: bool, by_cell: bool):

    assert rescale or by_cell

    scores['raw_score'] = scores.score

    if by_cell:
        pos_mean_scores_by_cell = scores[scores.score > 0].groupby(['func', 'cell_id']).score.mean().to_dict()
        neg_mean_scores_by_cell = scores[scores.score < 0].groupby(['func', 'cell_id']).score.mean().to_dict()
        scores['score'] = scores.apply(
            lambda r: (
                r.score / pos_mean_scores_by_cell[r.func, r.cell_id]
                if r.score > 0 else
                -r.score / neg_mean_scores_by_cell[r.func, r.cell_id]
            ),
            axis=1
        )

    if rescale:
        # feature scaling
        scores_by_func = scores.groupby('func').score
        scores_min = scores_by_func.min().to_dict()
        scores_range = (
            scores_by_func.max() - scores_by_func.min()
        ).to_dict()
        scores.score = scores.apply(lambda r: (
            -1 + (r.score - scores_min[r.func]) * 2 / scores_range[r.func]
        ), axis=1)

    return scores
