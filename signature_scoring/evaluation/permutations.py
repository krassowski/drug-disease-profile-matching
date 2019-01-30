import pickle
from copy import copy
from time import sleep
from types import FunctionType
from typing import List
from warnings import warn

from pandas import DataFrame, concat

from multiprocess import Pool

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

    for permutation, subtype in reevaluated_permutations:
        permutation['subtype'] = subtype

    return concat(
        permutation
        for permutation, subtype in reevaluated_permutations
    )
