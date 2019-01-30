import random
from collections import defaultdict
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO
from typing import Dict
from warnings import warn

from pandas import concat, DataFrame, Categorical, Series
from tqdm import tqdm

from .display import maximized_metrics, minimized_metrics, choose_columns


def subtypes_benchmark(expression, samples_by_subtype, benchmark_function, funcs, *args, samples_mapping=lambda x: x, **kwargs):
    subtypes_results = {}

    for type, samples in samples_by_subtype.items():
        type_subset = expression[[samples_mapping(sample) for sample in samples]]
        print(f'Using subset: {type} with {len(type_subset.columns)} samples')
        differential_subset = type_subset.differential('tumor', 'normal', only_paired=False)
        subtypes_results[type] = benchmark_function(funcs, differential_subset, *args, **kwargs)

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


def compare_against_permutations_group(
    function_result: Series, function_permutations: DataFrame,
    minimized_columns, maximized_columns,
    include_permutations=False

):
    data = []

    for metrics, sign in [(minimized_columns, -1), (maximized_columns, 1)]:

        for metric_name in metrics:
            observed_value = function_result[metric_name]

            if metric_name not in function_permutations.columns:
                warn(
                    f'Skipping {metric_name} (not present in permutations data). '
                    f'Please reevaluate permutations to include this metric'
                )
                continue

            metric_permutations = function_permutations[metric_name]

            more_extreme = (
                metric_permutations > observed_value
                if sign == 1 else
                metric_permutations < observed_value
            )
            p_value = sum(more_extreme) / len(metric_permutations)

            datum = {
                'p_value': p_value,
                'metric': metric_name,
                'observed': observed_value,
            }

            if include_permutations:
                datum['permutations'] = metric_permutations.tolist()

            data.append(datum)

    return DataFrame(data)


def compare_observations_with_permutations(
    subtypes_results,
    permutations,
    ranked_categories={'indications', 'contraindications', 'controls'},
    check_functions=True,
    reevaluate_permutations=False,
    reevaluation_kwargs=None
):

    permutations_grouped_by_corresponding_cluster = group_permutations_by_subtype(permutations)

    data = []

    for subtype, permutations in tqdm(permutations_grouped_by_corresponding_cluster.items()):
        result = subtypes_results[subtype]

        maximized_columns = choose_columns(result, maximized_metrics, ranked_categories)
        minimized_columns = choose_columns(result, minimized_metrics, ranked_categories)

        if check_functions:
            # do we have same scoring functions in permutations and observations?
            assert set(result.index) == set(permutations.index)

        for scoring_function in tqdm(permutations.index.unique()):
            function_result = result.loc[scoring_function]
            function_permutations = permutations.loc[scoring_function]

            rows = (
                compare_against_permutations_group(
                    function_result, function_permutations, minimized_columns, maximized_columns,
                    include_permutations=True, reevaluate_permutations=reevaluate_permutations,
                    reevaluation_kwargs=reevaluation_kwargs
                )
            )

            rows['subtype'] = subtype
            rows['scoring_function'] = scoring_function

            data.append(rows)

    df = concat(data)

    # lets save some memory
    categorical_variables = ['subtype', 'scoring_function', 'metric']

    # TODO: would be great to work on https://github.com/pandas-dev/pandas/issues/4464
    #  on some rainy weekend - this could have a huge speed & memory benefit
    for categorical_name in categorical_variables:
        df[categorical_name] = Categorical(df[categorical_name])

    return df

