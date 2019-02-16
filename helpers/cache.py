from pandas.util import hash_pandas_object
from pandas import DataFrame, Series
from copy import copy
from hashlib import sha512


def deep_hash(item):
    if isinstance(item, tuple):
        item = tuple(
            deep_hash(i)
            for i in item
        )
    if isinstance(item, dict):
        item = tuple(item.items())
    return hash(item)


def hash_series(series):
    """Fast series hash"""

    # hash(series.to_string())
    # 112 ms

    # hash(tuple(series.items()))
    # 2.24 ms

    # hash((tuple(series.values), tuple(series.index.values)))
    # 1.77 ms

    # chosen solution: 82.3 µs
    return hash((
        hash(series.index.values.data.tobytes()),
        hash(series.values.data.tobytes()),
    ))


def power_hash(obj):
    """Powerful (will hash almost everything)"""
    hashable = copy(obj)
    if isinstance(hashable, dict):
        hashable = hash(tuple(
            (key, power_hash(value))
            for key, value in obj.items()
        ))
    elif isinstance(hashable, DataFrame) or isinstance(hashable, Series):
        if isinstance(hashable, DataFrame):
            for column in hashable.columns:
                hashable[column] = hashable[column].apply(deep_hash)
        hashable = sha512(
            hash_pandas_object(hashable).values
        ).hexdigest()
    elif isinstance(hashable, list) or isinstance(hashable, tuple):
        hashable = tuple(
            power_hash(o)
            for o in obj
        )
    elif isinstance(hashable, set):
        hashable = frozenset(hashable)
    return hashable


verbose = 0


def cache_generator(make_hashable):

    def cache_decorator_closure(function):

        def cached(*args, **kwargs):
            hashable = make_hashable(args, kwargs)
            if not hasattr(function, '__cache__'):
                function.__cache__ = {}
            if hashable not in function.__cache__:
                if verbose > 1:
                    print(f'Adding to cache; key: {hashable}')
                function.__cache__[hashable] = function(*args, **kwargs)
            elif verbose > 1:
                print(f'Reusing cache; key: {hashable}')
            return function.__cache__[hashable]

        cached.original_function = function

        return cached

    return cache_decorator_closure


cache_decorator = cache_generator(lambda args, kwargs: power_hash((args, kwargs)))


def series_method_args_to_hashable(args, kwargs):
    return (
        hash_series(args[0]),   # series is assumed to be the first
        args[1:],
        tuple(kwargs.items())
    )


fast_series_cache_decorator = cache_generator(series_method_args_to_hashable)


def cached_property(function):
    
    def cached(self, *args, **kwargs):
        hashable = power_hash(self)
        if not hasattr(function, '__cache__'):
            function.__cache__ = {}
        if hashable not in function.__cache__:
            function.__cache__[hashable] = function(self, *args, **kwargs)
        return function.__cache__[hashable]
    
    return property(cached)
