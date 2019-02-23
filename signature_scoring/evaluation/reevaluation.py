from collections import defaultdict
from functools import partial
from types import SimpleNamespace
from warnings import warn

from numpy import isclose
from pandas import DataFrame, Series, Categorical
from pandas.core.dtypes.common import is_numeric_dtype

from data_frames import to_nested_dicts, the_only_one
from signature_scoring.evaluation import summarize_across_cell_lines, evaluation_summary
from signature_scoring.scoring_functions import ScoringFunction


def extract_single_score(g: DataFrame):
    return the_only_one(g.score)


def reevaluate(scores_dict_by_cell, scoring_func: ScoringFunction, subtypes_top=None, subtypes_scores=None, **kwargs):
    # TODO: support scoring functions with no cell grouping
    # if scoring_func.grouping:
        # pass
    summarize_test = partial(evaluation_summary, **kwargs)
    data = summarize_across_cell_lines(
        summarize_test, scores_dict_by_cell,
        subtypes_top=subtypes_top, subtypes_scores=subtypes_scores
    )
    return data


def extract_scores_from_result(
    result: Series,
    scores_as_series=True, are_grouped_by_cell=True
) -> DataFrame:

    if scores_as_series:
        def row_details(scores: Series) -> dict:
            yield {'score': scores}
    else:
        def row_details(scores: Series) -> dict:
            for row in scores.reset_index().itertuples():
                yield row._asdict()

    if are_grouped_by_cell:
        def iter_groups(scores: SimpleNamespace):
            for cell_id, cell_scores in scores.__dict__.items():
                yield {'cell_id': cell_id}, cell_scores
    else:
        def iter_groups(scores: SimpleNamespace):
            yield {}, scores

    data = []

    for func, func_scores in result.iteritems():
        for group_metadata, group_scores in iter_groups(func_scores):
            # group: indications / contra / controls
            for group, score_series in group_scores.__dict__.items():

                for row in row_details(score_series):
                    data.append({
                        'func': func,
                        'group': group,
                        **row,
                        **group_metadata
                    })

    return DataFrame(data)


def reevaluate_benchmark(old_benchmark: DataFrame, reevaluate_kwargs, verbose=True, keep_scores=True) -> DataFrame:
    # TODO: support scoring functions with no cell grouping
    scores = extract_scores_from_result(old_benchmark['meta:Scores'])

    if not any(scores_series.score.any() for scores_series in scores.score):
        warn(f'Skipping re-evaluation for {", ".join(old_benchmark.index)}: no scores in the original result')
        # return just empty df with same columns
        return DataFrame()

    scores_dict_by_func_cell_group_subtype = to_nested_dicts(
        scores, ['func', 'cell_id', 'group'],
        extract=extract_single_score
    )

    del scores

    data = []
    for func, scores_dict_by_cell_group in scores_dict_by_func_cell_group_subtype.items():
        reevaluate_kwargs['aggregate'] = None     # already aggregated
        result = reevaluate(
            scores_dict_by_cell_group,
            func,
            **reevaluate_kwargs
        )
        data.append({**result, **{'Func': func}})

    reevaluated_benchmark = DataFrame(data).set_index('Func')

    # preserve execution time information
    reevaluated_benchmark['Time'] = old_benchmark['Time']

    if verbose:

        old_columns = set(old_benchmark.columns)
        new_columns = set(reevaluated_benchmark.columns)

        removed_metrics = old_columns - new_columns
        if removed_metrics:
            print(removed_metrics, 'metrics removed')

        added_metrics = new_columns - old_columns
        if added_metrics:
            print(added_metrics, 'metrics added')

        old_benchmark_reordered = old_benchmark.loc[reevaluated_benchmark.index]
        for column in (old_columns & new_columns):
            if not old_benchmark_reordered[column].equals(reevaluated_benchmark[column]):
                if (
                    is_numeric_dtype(old_benchmark_reordered[column].dtype)
                    and
                    isclose(old_benchmark_reordered[column], reevaluated_benchmark[column]).all()
                ):
                    print(column, 'has changed slightly (float isclose)')
                else:
                    print(column, 'changed')

    if not keep_scores:
        reevaluated_benchmark.drop('meta:Scores', axis=1, inplace=True)

    return reevaluated_benchmark
