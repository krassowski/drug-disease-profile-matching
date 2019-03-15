from copy import copy
from glob import glob
from warnings import warn

from pandas import Categorical, concat, DataFrame

from signature_scoring.evaluation import permutations
from signature_scoring.evaluation.analysis import compute_metrics_by_func, normalize_scores
from signature_scoring.evaluation.reevaluation import extract_scores_from_result


def scores_from_subtypes(stratification_result, by_cell):
    # note: rescaling after merge! another option would be to rescale within subtypes,
    # giving each subtype equal deciding vote; but that would ignore the fact that they
    # may simply have little opinion on particular topic (i.e they do not know which drugs
    # would match and know about it) and force them to decide. On the other hand this approach
    # makes certain subtypes with greater number of cases more "powerful" for certain methods
    # which give "stronger" score when the number of analyzed cases is greater
    scores = concat([
        extract_scores_from_result(
            stratification_result[subtype]['meta:Scores'],
            scores_as_series=False,
            are_grouped_by_cell=by_cell
        ).assign(subtype=subtype)
        for subtype in stratification_result
    ]).reset_index(drop=True)
    scores = normalize_scores(scores, rescale=True, by_cell=by_cell)
    scores['is_indication'] = scores.group == 'indications'
    scores['is_indication'] = Categorical(scores.is_indication)
    return scores


def get_all_scores(results, by_cell):
    scores_all_stratifications = []
    for stratification, stratification_result in results.items():

        results_with_scores = {}
        for subtype in stratification_result:
            if 'meta:Scores' in stratification_result[subtype]:
                results_with_scores[subtype] = stratification_result[subtype]
            else:
                warn(f'{subtype} of {stratification} has no scores for this permutation')

        scores = scores_from_subtypes(results_with_scores, by_cell=by_cell)
        scores_all_stratifications.append(
            scores.assign(stratification=stratification)
        )
    return scores_all_stratifications


def aggregate_scores_per_substance_strat_func(scores_with_reference, trans='max'):
    by = ['func', 'pert_iname', 'stratification']

    df: DataFrame = copy(scores_with_reference)

    transformed_scores = df.groupby(by).score.transform(trans)
    df.score = transformed_scores

    df = df[df.columns.difference(['pert_idose', 'cell_id', 'subtype', 'raw_score'])]
    df = df.drop_duplicates()

    return df


def metrics_by_stratification(scores_with_stratification):
    pan_stratification_metrics = []
    for stratification, scores in scores_with_stratification.groupby('stratification'):
        metrics = compute_metrics_by_func(scores, without_unassigned=True)
        pan_stratification_metrics.append(
            metrics.assign(stratification=stratification)
        )
    return concat(pan_stratification_metrics)


def rename_for_plot(data, function_names, stratification_names):
    data.func = data.func.replace(function_names)
    data.func = Categorical(
        data.func,
        ordered=True,
        categories=function_names.values()
    )
    data.stratification = data.stratification.replace(stratification_names)
    data.stratification = Categorical(
        data.stratification,
        ordered=True,
        categories=stratification_names.values()
    )
    return data


def prepare_metrics_for_plotting(metrics, nice_names, stratification_nice_names):
    stratification_metrics = metrics.melt(id_vars=['func', 'stratification'])
    return rename_for_plot(
        stratification_metrics, nice_names, stratification_nice_names
    )


def load_pickled_permutations(path='../../data/subtype_permutations/', prefix=''):
    all_permutations = []
    for file in glob(path + prefix + '*.pickle'):
        try:
            all_permutations.extend(permutations.load(file.replace('.pickle', '')))
        except:
            warn(f'{file} permutations broken')
    return all_permutations
