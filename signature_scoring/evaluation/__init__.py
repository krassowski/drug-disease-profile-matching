from types import FunctionType
from typing import Dict

import numpy
from functools import partial

from pandas import DataFrame
from operator import attrgetter

from data_sources.drug_connectivity_map import Scores, dcm
from helpers.gui.namespace import NeatNamespace
from helpers.cache import cache_decorator
from helpers import WarningManager

from .. import score_signatures
from ..models import SignaturesGrouping
from .scores_models import ScoresVector, ProcessedScores, TopScores
from .metrics import EvaluationMetric, metrics_manager

test_warnings = WarningManager()


def evaluation_summary(scores_dict, top, aggregate, ScoresVector=ScoresVector, transform_scores=None):
    if aggregate:
        aggregate = attrgetter(aggregate)
        scores_dict = {label: aggregate(scores) for label, scores in scores_dict.items()}

    if transform_scores:
        scores_dict = transform_scores(scores_dict)

    scores_indications, scores_controls, scores_contraindications = (
        scores_dict['indications'], scores_dict['controls'], scores_dict['contraindications']
    )

    scores_map = {
        1: scores_indications,
        0: scores_controls,
        -1: scores_contraindications
    }

    vector = ScoresVector(scores_map)

    if top == 'rescaled':
        test_warnings.warn_once('Using 0.5 threshold of scaled scores vector to select top results.')
        top_scoring = vector.observed[vector.observed > 0.5]
    else:
        test_warnings.warn_once('Using top 0.1 quantile to select top results.')
        top_scoring = vector.observed[vector.observed > vector.observed.quantile(0.9)]

    scores = ProcessedScores(
        vector_overall=vector,
        vector_contraindications=ScoresVector(scores_map, limit_to=[1, -1]),
        vector_controls=ScoresVector(scores_map, limit_to=[1, 0]),
        top=TopScores(
            all=top_scoring,
            indications=top_scoring[scores_indications.index.intersection(top_scoring.index)],
            contraindications=top_scoring[scores_contraindications.index.intersection(top_scoring.index)],
            controls=top_scoring[scores_controls.index.intersection(top_scoring.index)]
        ),
        **{
            category: aggregated_scores.score
            for category, aggregated_scores in scores_dict.items()
        }
    )

    results = {
        'meta:Selected substances': len(top_scoring),
        'meta:Scores': NeatNamespace(scores_dict),
    }

    for category, metrics in metrics_manager.registry.items():
        category_scores = scores_dict.get(category, [])
        if category == 'overall' or len(category_scores):
            for metric in metrics.values():
                key = f'{category}:{metric.name}'
                assert key not in results
                results[key] = metric(scores)

    return results


#@cache_decorator
def calculate_scores(query_signature, signatures_map, scoring_func, fold_changes, **kwargs):
    score = partial(
        score_signatures,
        scoring_func, query_signature,
        fold_changes=fold_changes, warning_manager=test_warnings, **kwargs
    )

    def score_or_empty(category):
        if category not in signatures_map:
            return Scores({})
        signatures = signatures_map[category]
        return score(signatures)

    scores_dict = {
        'indications': score(signatures_map['indications']),
        'controls': score_or_empty('control'),
        'contraindications': score_or_empty('contraindications')
    }

    return scores_dict


def select_cells(signatures_map: Dict[str, SignaturesGrouping], completeness_ratio: float):
    all_signatures = {
        signature
        for grouping in signatures_map.values()
        for signature in grouping.signature_ids
    }
    all_data = dcm.sig_info[dcm.sig_info.sig_id.isin(all_signatures)]
    count_by_cell = all_data.drop_duplicates(['cell_id', 'pert_iname']).groupby('cell_id').count().pert_iname

    all_substances_and_controls_cnt = len(set(all_data.pert_iname))

    selected_cells_counts = count_by_cell[count_by_cell >= all_substances_and_controls_cnt * completeness_ratio]
    selected_cells = set(selected_cells_counts.index)

    selected_cells = {
        cell
        for cell in selected_cells
        if all(
            cell in set(all_data[all_data.sig_id.isin(grouping.signature_ids)].cell_id)
            for grouping in signatures_map.values()
        )
    }
    return selected_cells


def combine_values(column):
    by_categories = metrics_manager.registry
    category, metric_id = column.name.split(':')
    try:
        metric = by_categories[category][metric_id]
        return metric.combine(column.tolist())
    except KeyError:
        try:
            if category != 'meta':
                print(
                    f'Failed combining to determine preferred way to combine'
                    f'values in column: {column.name}, trying mean() instead'
                )
            return column.mean()
        except TypeError:
            return numpy.nan


def summarize_across_cell_lines(summarize_test: FunctionType, scores_dict_by_cell: dict):
    data = []
    for cell_id, cell_scores_dict in scores_dict_by_cell.items():
        summary = summarize_test(cell_scores_dict)
        summary['meta:cell_id'] = cell_id
        data.append(summary)
    data = DataFrame(data)
    scores = NeatNamespace(data.set_index('meta:cell_id')['meta:Scores'].to_dict())

    data = data.apply(combine_values)
    data = data.to_dict()
    data['meta:Scores'] = scores
    return data


def evaluate(
    scoring_func, query_signature, indications_signatures, contraindications_signatures,
    control_signatures=None, aggregate='mean_per_substance_dose_and_cell', top='rescaled',
    cell_lines_ratio=0.9, summary='per_cell_line_combined', fold_changes=False,
    reset_warnings=True, **kwargs
):
    """
    aggregate: mean_per_substance, best_per_substance, signal_to_noise
    """
    if reset_warnings:
        test_warnings.reset()

    signatures_map = {
        'indications': indications_signatures,
        'contraindications': contraindications_signatures,
        'control': control_signatures
    }

    if control_signatures is not None:
        if fold_changes:
            print('Skipping controls as those will be used for fold_change calculation')
            del signatures_map['control']

    collection = scoring_func.collection

    # remove signatures which were not given
    signatures_map: Dict[str, SignaturesGrouping] = {
        label: collection(signatures)
        for label, signatures in signatures_map.items()
        if signatures is not None
    }

    if cell_lines_ratio:

        selected_cells = select_cells(signatures_map, cell_lines_ratio)

        test_warnings.warn_once(
            f'Keeping cell lines with data for at least {cell_lines_ratio * 100:.0f}% '
            f'of substances (and at least one substance per class).'
        )
        test_warnings.warn_once(
            f'Keeping {len(selected_cells)} distinct cell lines: {selected_cells}'
        )

        for name, signatures in signatures_map.items():
            signatures_data = dcm.sig_info[dcm.sig_info.sig_id.isin(signatures.signature_ids)]
            signatures_to_keep = set(signatures_data[signatures_data.cell_id.isin(selected_cells)].sig_id)
            signatures_map[name] = signatures.drop_signatures(
                [column for column in signatures.signature_ids if column not in signatures_to_keep]
            )

    scores_dict = calculate_scores(query_signature, signatures_map, scoring_func, fold_changes, **kwargs)

    summarize_test = partial(evaluation_summary, top=top, aggregate=aggregate)

    if cell_lines_ratio and summary == 'per_cell_line_combined' and not scoring_func.grouping:
        scores_dict_by_cell = {}

        for cell_id in selected_cells:
            cell_scores_dict = {
                label: scores.limit_to_cell_line(cell_id)
                for label, scores in scores_dict.items()
            }
            scores_dict_by_cell[cell_id] = cell_scores_dict

        data = summarize_across_cell_lines(summarize_test, scores_dict_by_cell)
    else:
        data = summarize_test(scores_dict)

    return data
