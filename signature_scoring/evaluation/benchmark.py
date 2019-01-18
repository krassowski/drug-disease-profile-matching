import gc
from time import time

from pandas import DataFrame
from tqdm import tqdm_notebook
from ..models import ExpressionWithControls
from . import evaluate


def benchmark(
    funcs, query_signature, indications_signatures, contraindications_signatures=None,
    control_signatures=None, per_test_progress=False, query_expression: ExpressionWithControls = None,
    quiet=False, progress=True, **kwargs
):
    data = []
    is_first_run = True
    if progress:
        funcs = tqdm_notebook(funcs)
    for func in funcs:
        if not quiet:
            print(f'Testing {func.__name__}')

        query = query_expression if func.input == ExpressionWithControls else query_signature

        start = time()
        result = evaluate(
            func, query, indications_signatures, contraindications_signatures,
            control_signatures=control_signatures if func.is_applicable_to_control_signatures else None,
            progress=per_test_progress, reset_warnings=is_first_run,
            **kwargs
        )
        end = time()
        data.append({**result, **{'Func': func.__name__, 'Time': end - start}})

        gc.collect()
        is_first_run = False

    return DataFrame(data).set_index('Func')
