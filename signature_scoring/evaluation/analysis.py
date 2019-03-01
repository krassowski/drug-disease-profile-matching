from pandas import DataFrame, Series
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import log_loss

from helpers.gui import NeatNamespace
from helpers.r import r_ks_test

from .calculation_utilities import generalized_roc_auc_score


def ks_distance(a, b):
    return r_ks_test(a, b)['statistic']


def distance_matrix(scores, normalize=True, distance=ks_distance, condensed=True):
    functions = list(scores.func.unique())
    matrix = DataFrame(data=0, index=functions, columns=functions)
    for i, func in enumerate(functions):
        for other_func in functions[i+1:]:
            a = scores[scores.func == func]['score']
            b = scores[scores.func == other_func]['score']
            d = distance(a, b)
            matrix.loc[func, other_func] = d
            if not condensed:
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

    func_masks = {
        func: scores.func == func
        for func in scores.func.unique()
    }

    if by_cell:
        pos_mean_scores_by_cell = scores[scores.score > 0].groupby(['func', 'cell_id']).score.mean().to_dict()
        neg_mean_scores_by_cell = scores[scores.score < 0].groupby(['func', 'cell_id']).score.mean().to_dict()

        # 2.49 s
        positive = scores.score > 0
        s_loc = 'score'

        cell_masks = {
            cell_id: scores.cell_id == cell_id
            for cell_id in scores.cell_id.unique()
        }

        for sign, reference, leg in [
            (+1, pos_mean_scores_by_cell, positive),
            (-1, neg_mean_scores_by_cell, ~positive)
        ]:
            for func, mask in func_masks.items():
                leg_and_func = leg & mask
                for cell_id in scores.cell_id.unique():
                    denominator = reference[func, cell_id] * sign
                    loc = leg_and_func & cell_masks[cell_id]
                    scores.loc[loc, s_loc] = scores.loc[loc, s_loc] / denominator

        # 4.63 s
        # scores['score'] = scores.apply(
        #     lambda r: (
        #         r.score / pos_mean_scores_by_cell[r.func, r.cell_id]
        #         if r.score > 0 else
        #         -r.score / neg_mean_scores_by_cell[r.func, r.cell_id]
        #     ),
        #     axis=1
        # )

    if rescale:
        # feature scaling
        scores_by_func = scores.groupby('func').score
        scores_min = scores_by_func.min().to_dict()
        scores_range = (
            scores_by_func.max() - scores_by_func.min()
        ).to_dict()

        # 303 ms
        for func, mask in func_masks.items():
            scores.loc[mask, 'score'] = (
                -1 + (scores.loc[mask, 'score'] - scores_min[func]) * 2 / scores_range[func]
            )

        # 2.49 s
        # scores.score = scores.apply(lambda r: (
        #    -1 + (r.score - scores_min[r.func]) * 2 / scores_range[r.func]
        # ), axis=1)

    return scores


def test_rank_diff_indications_vs_non_indications(scores, alternative_indications_are='greater'):
    assert 'rank' in scores.columns
    is_indication = (scores.group == 'indications')
    indications_by_func = scores[is_indication].groupby('func')
    non_indications_by_func = scores[~is_indication].groupby('func')
    non_indications_by_func = non_indications_by_func['rank'].apply(list).to_dict()
    result = indications_by_func['rank'].apply(
        lambda ranks: r_ks_test(
            ranks,
            Series(non_indications_by_func[ranks.name]),
            alternative=alternative_indications_are
        )
    )
    result = result.reset_index().set_index('func').pivot(columns='level_1')
    result.columns = result.columns.droplevel()
    result.drop('data.name', axis=1, inplace=True)
    return result


def n_best_scoring(scores, group_by_function, n=3):
    best_scoring_pos = (
        group_by_function['rank']
        .apply(Series.nsmallest, n)
        .reset_index()
        .rename({'level_1': 'score_index'}, axis=1)
    )
    best_recovered = scores.loc[best_scoring_pos.score_index].set_index('func')
    return best_recovered


def clusters_from_matrix(matrix, n, linkage=None):
    """If no specific linkage is given, compare results
    from all implemented ones and check if they agree"""
    if linkage:
        linkages = [linkage]
    else:
        linkages = ["ward", "complete", "average", "single"]

    clusterings = []
    for linkage in linkages:
        clusters = list(
            AgglomerativeClustering(n_clusters=n, linkage=linkage)
            .fit_predict(matrix)
        )
        clusterings.append(clusters)

    try:
        assert all(adjusted_rand_score(clusters, c) for c in clusterings)
    except:
        print(clusterings)

    clusters = dict(zip(matrix.columns, clusters))
    return clusters
