from types import FunctionType
from typing import Dict, Union

import numpy
from functools import partial

from pandas import DataFrame, concat, Series
import pandas
from operator import attrgetter

from data_frames import is_copy
from data_sources.drug_connectivity_map import Scores, dcm, AggregatedScores
from helpers.gui.namespace import NeatNamespace
from helpers import WarningManager

from .. import score_signatures
from ..models import SignaturesGrouping
from .scores_models import ScoresVector, ProcessedScores, TopScores, Group
from .metrics import EvaluationMetric, metrics_manager


pandas.options.mode.chained_assignment = None

test_warnings = WarningManager()


groups_label_value_map = {
    'indications': 1,
    'controls': 0,
    'contraindications': -1
}


def transform_to_aggregated_scores_map(
    scores_by_group: Dict[Group, Scores], aggregate: str,
    transform_scores=None
) -> Dict[str, AggregatedScores]:

    def verify_aggregated(scores):
        assert isinstance(scores, AggregatedScores)
        return scores

    aggregate = (
        attrgetter(aggregate)
        if aggregate else
        verify_aggregated
    )

    scores_by_group = {label: aggregate(scores) for label, scores in scores_by_group.items()}

    if transform_scores:
        scores_by_group = transform_scores(scores_by_group)

    return {
        group: scores
        for group, scores in scores_by_group.items()
    }


def select_top_substance(vector: ScoresVector, how) -> Series:
    if how == 'rescaled':
        test_warnings.warn_once('Using 0.5 threshold of scaled scores vector to select top results.')
        top_scoring = vector.observed[vector.observed > 0.5]
    elif how == 'quantile':
        test_warnings.warn_once('Using top > 0.1 quantile to select top results.')
        top_scoring = vector.observed[vector.observed > vector.observed.quantile(0.9)]
    else:
        assert False
    return top_scoring


def evaluation_summary(
    scores_dict: Dict[Group, Union[Scores, AggregatedScores]], top: str, aggregate: str, scores_vector=ScoresVector,
    transform_scores=None, subtypes_top=None, subtypes_dicts=None
):
    # TODO: rename top -> top_selection_criteria

    aggregated_scores_by_group_label = transform_to_aggregated_scores_map(scores_dict, aggregate, transform_scores)

    aggregated_scores_df = concat(
        aggregated_scores_by_group_label[group].assign(group=group, expected_score=value)
        for group, value in groups_label_value_map.items()
    )
    aggregated_scores_df['group'] = pandas.Categorical(aggregated_scores_df['group'])

    vector = scores_vector(aggregated_scores_df)

    if subtypes_top:

        assert subtypes_top == 'best_from_each_subtype'

        top_scoring_in_subtypes = []

        for subtype, subtype_scores_dict in subtypes_dicts.items():
            subtype_map = transform_to_aggregated_scores_map(subtype_scores_dict, aggregate, transform_scores)

            # contaminate data frames by in-place assignment to reduce copy operations
            for group, value in groups_label_value_map.items():
                # this is a dict selection, no need to check that copy gets returned
                # (plus I assume that it will change the data - with little impact)
                df = subtype_map[group]
                df['group'] = group
                df['expected_score'] = value

            aggregated_subtype_scores = concat(subtype_map.values())
            subtype_vector = scores_vector(aggregated_subtype_scores)
            top_scoring_in_subtypes.append(
                select_top_substance(subtype_vector, how=top)
            )

        top_scoring = concat(top_scoring_in_subtypes)
        top_scoring = top_scoring.groupby(top_scoring.index.names).max()

    else:
        top_scoring = select_top_substance(vector, how=top)

    scores_indications, scores_controls, scores_contraindications = (
        aggregated_scores_by_group_label['indications'],
        aggregated_scores_by_group_label['controls'],
        aggregated_scores_by_group_label['contraindications']
    )

    top = TopScores(
        all=top_scoring,
        indications=top_scoring[scores_indications.index.intersection(top_scoring.index)],
        contraindications=top_scoring[scores_contraindications.index.intersection(top_scoring.index)],
        controls=top_scoring[scores_controls.index.intersection(top_scoring.index)]
    )

    top_set = set(top_scoring.index.to_series())

    def one_if_in_top(substance):
        return 1 if substance in top_set else 0

    aggregated_with_score_by_top = aggregated_scores_df.assign(
        score=Series([one_if_in_top(substance) for substance in aggregated_scores_df.index]).values
    )

    def rescore(expected, df=aggregated_with_score_by_top):
        df_subset = df[df.group.isin(expected)]
        assert is_copy(df_subset, df)
        df_subset['expected_score'] = df_subset.group.map(expected.get)
        return df_subset

    indications_over_controls = rescore(expected={'indications': 1, 'controls': 0})
    indications_over_contraindications = rescore(expected={'indications': 1, 'contraindications': 0})
    indications_over_controls_and_contra = rescore(expected={'indications': 1, 'contraindications': 0, 'controls': 0})

    scores = ProcessedScores(
        vector_overall=vector,
        vector_contraindications=scores_vector(aggregated_scores_df, limit_to=['indications', 'contraindications'], rescale=False),
        vector_controls=scores_vector(aggregated_scores_df, limit_to=['indications', 'controls'], rescale=False),
        vector_overall_binary=scores_vector(indications_over_controls_and_contra, rescale=False),
        vector_contraindications_binary=scores_vector(indications_over_contraindications, rescale=False),
        vector_controls_binary=scores_vector(indications_over_controls, rescale=False),
        top=top,
        **{
            group: aggregated_scores.score
            for group, aggregated_scores in aggregated_scores_by_group_label.items()
        }
    )

    results = {
        'meta:Selected substances': len(top_scoring),
        'meta:Scores': NeatNamespace(aggregated_scores_by_group_label),
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


def summarize_across_cell_lines(summarize_test: FunctionType, scores_dict_by_cell: dict, subtypes_top=None, subtypes_scores=None):
    data = []
    for cell_id, cell_scores_dict in scores_dict_by_cell.items():
        subtype_score = None if not subtypes_scores else subtypes_scores[cell_id]
        summary = summarize_test(cell_scores_dict, subtypes_top=subtypes_top, subtypes_dicts=subtype_score)
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

    # remove signatures that were not given
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
