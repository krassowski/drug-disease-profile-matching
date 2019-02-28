import random
from collections import defaultdict
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO
from typing import Dict
from warnings import warn

from pandas import concat, DataFrame, Categorical
from tqdm import tqdm

from data_frames import to_nested_dicts
from data_sources.drug_connectivity_map import AggregatedScores
from ..models.with_controls import TCGAExpressionWithControls

from .permutations import compare_against_permutations_group, compare_observations_with_permutations
from .reevaluation import extract_scores_from_result, extract_single_score, reevaluate
from .display import maximized_metrics, minimized_metrics, choose_columns
from .scores_models import Group


def subtypes_benchmark(
    expression, samples_by_subtype, benchmark_function, funcs, *args,
    samples_mapping=lambda x: x, use_all_controls=True,
    single_sample=True, multi_sample=True, **kwargs
):
    subtypes_results = {}

    if use_all_controls:
        all_controls = expression[expression.columns[expression.classes == 'normal']]
        additional_controls = {'additional_controls': all_controls}
    else:
        additional_controls = {}

    print(f'Using all disease controls: {use_all_controls}')

    for subtype, samples in samples_by_subtype.items():

        type_subset = expression[[samples_mapping(sample) for sample in samples]]

        queries = {'query_signature': None}

        if single_sample:
            print(f'Using subset: {subtype} with {len(type_subset.columns)} samples')
            if 'normal' not in type_subset.classes and not use_all_controls:
                print(
                    'No normal subtype-specific samples to create differential expression, '
                    'set use_all_controls=True to include control samples from all subtypes.'
                )
                continue
            differential_subset = type_subset.differential(
                'tumor', 'normal',
                only_paired=False,
                **additional_controls
            )
            if differential_subset is None:
                print(f'Skipping subtype {subtype}')
                continue
            queries['query_signature'] = differential_subset

        if multi_sample:
            if use_all_controls:
                absent_controls = all_controls.columns.difference(type_subset.columns)
                type_subset = concat([type_subset, all_controls[absent_controls]], axis=1)

            subset_with_controls = TCGAExpressionWithControls(type_subset)
            queries['query_expression'] = subset_with_controls

        subtypes_results[subtype] = benchmark_function(
            funcs,
            *args, **{**queries, **kwargs}
        )

    return subtypes_results


def random_subtypes_benchmark(i, expression, *args, **kwargs):
    samples = list(expression.columns)
    random_mapping = dict(zip(samples, random.sample(samples, len(samples))))
    f = StringIO()
    with redirect_stdout(f), redirect_stderr(f):
        result = subtypes_benchmark(expression, *args, samples_mapping=random_mapping.get, **kwargs)
    return result


def group_permutations_by_subtype(permutations) -> Dict[str, DataFrame]:
    grouped_by_corresponding_cluster = defaultdict(list)

    for permutation in permutations:
        for cluster_name, result in permutation.items():
            grouped_by_corresponding_cluster[cluster_name].append(result)

    grouped_by_corresponding_cluster = {
        cluster_name: concat(results)
        for cluster_name, results in grouped_by_corresponding_cluster.items()
    }

    return grouped_by_corresponding_cluster


def test_permutations_number_in_subtype(
    result_subtype_subset: DataFrame, subtype_subset: DataFrame, subtype: str,
    ranked_categories={'indications', 'contraindications', 'controls'}
):
    data = []

    maximized_columns = choose_columns(result_subtype_subset, maximized_metrics, ranked_categories)
    minimized_columns = choose_columns(result_subtype_subset, minimized_metrics, ranked_categories)

    if set(subtype_subset.index.unique()) != set(result_subtype_subset.index.unique()):
        warn('Different sets of functions in result and permutations')

    for scoring_function in subtype_subset.index.unique():

        function_subset = subtype_subset.loc[scoring_function]
        result_function_subset = result_subtype_subset.loc[scoring_function]
        count = len(function_subset)

        for first_n_permutations in range(2, count + 1):
            considered_permutations = function_subset.head(first_n_permutations)

            measurements = DataFrame(compare_against_permutations_group(
                result_function_subset, considered_permutations,
                minimized_columns, maximized_columns
            ))
            measurements['n_permutations'] = first_n_permutations
            measurements['scoring_function'] = scoring_function
            data.append(measurements)

    joined = concat(data)
    joined['subtype'] = subtype
    return joined


def compare_observations_with_subtype_permutations(
    subtypes_results: Dict[Group, DataFrame],
    permutations: DataFrame,
    ranked_categories={'indications', 'contraindications', 'controls'},
    check_functions=True
):

    permutations_grouped_by_corresponding_cluster = permutations.groupby('subtype')

    data = []

    for subtype, permutations in tqdm(permutations_grouped_by_corresponding_cluster):
        result = subtypes_results[subtype]

        function_results = compare_observations_with_permutations(
            result, permutations,
            ranked_categories, check_functions
        )

        for function_result in function_results:
            function_result['subtype'] = subtype
            data.append(function_result)

    df = concat(data)

    # lets save some memory
    categorical_variables = ['subtype', 'scoring_function', 'metric']

    # TODO: would be great to work on https://github.com/pandas-dev/pandas/issues/4464
    #  on some rainy weekend - this could have a huge speed & memory benefit
    for categorical_name in categorical_variables:
        df[categorical_name] = Categorical(df[categorical_name])

    return df


def merge_subtypes_results_and_reevaluate(result, strategy, top):
    assert strategy in {'equal_weight_for_subtypes', 'weight_by_actual_score'}

    should_give_equal_weight_to_subtypes = (strategy == 'equal_weight_for_subtypes')

    extracted_scores = []
    for subtype, subtype_result in result.items():
        df = extract_scores_from_result(subtype_result['meta:Scores'])
        df['subtype'] = subtype
        extracted_scores.append(df)
    extracted_scores = concat(extracted_scores)

    scores_dict_by_func_cell_group_subtype = to_nested_dicts(extracted_scores, ['func', 'cell_id', 'group'])

    data = []

    for func, scores_dict_by_cell_group_subtype in scores_dict_by_func_cell_group_subtype.items():
        scores_dict_by_cell = {}
        scores_dict_by_cell_subtype_group = to_nested_dicts(
            extracted_scores.query('func == @func'),
            ['cell_id', 'subtype', 'group'],
            extract=extract_single_score
            # dict_type=NeatNamespace - for visual debugging
        )

        for cell_id, cell_scores in scores_dict_by_cell_group_subtype.items():
            cell_dict = {}
            for group, scores_by_subtype in cell_scores.items():
                score_series = concat(scores_by_subtype.score.tolist())
                index = score_series.index
                # take the best score by substance.
                cell_dict[group] = AggregatedScores(score_series.groupby(level=list(range(index.nlevels))).max())

            scores_dict_by_cell[cell_id] = cell_dict

        result = reevaluate(
            scores_dict_by_cell,
            func,
            top=top,
            aggregate=None,  # already aggregated
            subtypes_scores=scores_dict_by_cell_subtype_group,
            subtypes_top=(
                'best_from_each_subtype'
                if should_give_equal_weight_to_subtypes else
                None  # weight_by_actual_score is the default action
            )
        )
        data.append({**result, **{'Func': func}})

    return DataFrame(data).set_index('Func')
