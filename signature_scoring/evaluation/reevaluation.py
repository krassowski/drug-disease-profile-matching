from functools import partial
from warnings import warn

from numpy import isclose
from pandas import DataFrame
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


def extract_scores_from_result(result: DataFrame) -> DataFrame:
    # TODO: support scoring functions with no cell grouping
    data = []
    for func, scores in result['meta:Scores'].iteritems():
        for cell_id, cell_scores in scores.__dict__.items():
            # group: indications / contra / controls
            for group, score_series in cell_scores.__dict__.items():
                data.append({
                    'func': func,
                    'cell_id': cell_id,
                    'group': group,
                    'score': score_series
                })
    return DataFrame(data)


def reevaluate_benchmark(old_benchmark: DataFrame, reevaluate_kwargs, verbose=True, keep_scores=True) -> DataFrame:
    # TODO: support scoring functions with no cell grouping
    scores = extract_scores_from_result(old_benchmark)

    if not any(scores_series.score.any() for scores_series in scores.score):
        warn(f'Skipping re-evaluation for {", ".join(old_benchmark.index)}: no scores in the original result')
        # return just empty df with same columns
        return DataFrame()

    scores_dict_by_func_cell_group_subtype = to_nested_dicts(
        scores, ['func', 'cell_id', 'group'],
        extract=extract_single_score
    )

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

        for column in (old_columns & new_columns):
            if not old_benchmark[column].equals(reevaluated_benchmark[column]):
                if (
                    is_numeric_dtype(old_benchmark[column].dtype)
                    and
                    isclose(old_benchmark[column], reevaluated_benchmark[column]).all()
                ):
                    print(column, 'has changed slightly (float isclose)')
                else:
                    print(column, 'changed')

    if not keep_scores:
        reevaluated_benchmark.drop('meta:Scores', axis=1, inplace=True)

    return reevaluated_benchmark
