import pickle
from copy import copy
from time import sleep
from types import FunctionType
from typing import List
from warnings import warn

from pandas import DataFrame, concat, Series
from tqdm.auto import tqdm


from multiprocess import Pool

from .display import choose_columns, maximized_metrics, minimized_metrics
from .reevaluation import reevaluate_benchmark


def generate(
    randomizer: FunctionType, expression, samples_by_type, benchmark_partial, funcs,
    *signatures, n=100, pickle_name=None, processes=None,
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
    permutations = list(pool.imap(randomizer, range(n), args))

    if pickle_name:
        with open(f'{pickle_name}.pickle', 'wb') as f:
            pickle.dump(permutations, f)

    return permutations


def load(name):
    with open(f'{name}.pickle', 'rb') as f:
        permutations = pickle.load(f)
        return permutations


def ensure_kwargs(kwargs):

    if not kwargs:
        warn(
            'No kwargs for evaluate() provided - using the default top selection strategy '
            '(and other benchmarking parameters)'
        )
    assert kwargs is not None, 'You must provide revaluation kwargs (like top)'


def reevaluate(permutations: DataFrame, processes=None, **kwargs):
    """Some permutations were evaluated when not all the evaluation metrics were defined,

    so those need re-evaluation to include missing metric's values"""

    ensure_kwargs(kwargs)

    # reevaluate rows separately, as Func values are not-unique (by permutation definition)
    reevaluated_permutations = Pool(processes).imap(
        reevaluate_benchmark,
        [permutations.iloc[[i]] for i in range(len(permutations))],
        shared_args=(
            kwargs,  # reevaluate kwargs
            False    # verbose=False
        )
    )

    return concat(reevaluated_permutations)


def pass_metadata_through(func):
    def wrapped(joined_chunk, *args, **kwargs):
        real_chunk, metadata = joined_chunk
        return func(real_chunk, *args, **kwargs), metadata
    return wrapped


def reevaluate_with_subtypes(permutations: List[DataFrame], processes=None, **kwargs):
    """Like reevaluate() but accepting unprocessed list of permutations,
    with the subtypes information not yet assigned.

    This function was implemented as as special case
    due to high memory usage of reevaluate()
    """

    ensure_kwargs(kwargs)

    reevaluated_permutations = Pool(processes).imap(
        pass_metadata_through(reevaluate_benchmark),
        [
            (result, subtype)
            for result_by_subtype in permutations
            for subtype, result in result_by_subtype.items()
        ],
        shared_args=(
            kwargs,  # reevaluate kwargs
            False    # verbose=False
        )
    )

    return concat(
        permutation.assign(subtype=subtype)
        for permutation, subtype in reevaluated_permutations
        if not permutation.empty
    )


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
    result: DataFrame,
    permutations: DataFrame,
    ranked_categories={'indications', 'contraindications', 'controls'},
    check_functions=True
):
    """result = reference result, observed result"""
    data = []

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
                include_permutations=True
            )
        )

        rows['scoring_function'] = scoring_function

        data.append(rows)

    return data
