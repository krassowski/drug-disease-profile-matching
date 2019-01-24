import pickle
import random
from collections import defaultdict
from contextlib import redirect_stdout, redirect_stderr
from copy import copy
from io import StringIO
from time import sleep
from typing import Dict

from pandas import concat, DataFrame, Categorical
from tqdm import tqdm

from multiprocess import Pool
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


def generate_permutations(
    expression, samples_by_type, benchmark_partial, funcs,
    *signatures, n=100, pickle_name=None, processes=None
):
    print('To abort permutations generation start, send keyboard interrupt now - waiting 5s')
    sleep(5)

    pool = Pool(processes)
    single_process_benchmark = copy(benchmark_partial)
    single_process_benchmark.keywords['processes'] = 1
    single_process_benchmark.keywords['progress'] = False

    args = [
        expression,
        samples_by_type, single_process_benchmark,
        funcs,
        *signatures
    ]
    permutations = list(pool.imap(random_subtypes_benchmark, range(n), args))

    if pickle_name:
        with open(f'{pickle_name}.pickle', 'wb') as f:
            pickle.dump(permutations, f)

    return permutations


def load_permutations(name):
    with open(f'{name}.pickle', 'rb') as f:
        permutations = pickle.load(f)
        return permutations


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
    function_result, function_permutations,
    minimized_columns, maximized_columns,
    include_permutations=False,
    **kwargs
):
    data = []

    for metrics, sign in [(minimized_columns, -1), (maximized_columns, 1)]:

        for metric_name in metrics:
            observed_value = function_result[metric_name]
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
                **kwargs
            }

            if include_permutations:
                datum['permutations'] = metric_permutations.tolist()

            data.append(datum)

    return data


def compare_observations_with_permutations(
    subtypes_results,
    permutations,
    ranked_categories={'indications', 'contraindications', 'controls'},
    check_functions=True
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

        for scoring_function in permutations.index.unique():
            function_result = result.loc[scoring_function]
            function_permutations = permutations.loc[scoring_function]

            rows = (
                compare_against_permutations_group(
                    function_result, function_permutations, minimized_columns, maximized_columns,
                    include_permutations=True,
                    subtype=subtype, scoring_function=scoring_function
                )
            )

            data.extend(rows)

    df = DataFrame(data)

    # lets save some memory
    categorical_variables = ['subtype', 'scoring_function', 'metric']

    # TODO: would be great to work on https://github.com/pandas-dev/pandas/issues/4464
    #  on some rainy weekend - this could have a huge speed & memory benefit
    for categorical_name in categorical_variables:
        df[categorical_name] = Categorical(df[categorical_name])

    return df

