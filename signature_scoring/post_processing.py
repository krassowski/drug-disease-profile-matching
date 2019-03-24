import gc

from pandas import Categorical

from .evaluation.analysis import normalize_scores


def rank_by_func(scores):
    scores['rank'] = scores.groupby('func').score.rank(ascending=False, method='first')
    return scores


def add_helper_columns(scores):
    scores['is_indication'] = scores.group == 'indications'
    scores['is_indication'] = Categorical(scores.is_indication)
    scores['known_indication'] = scores['is_indication'].cat.rename_categories({
        True: 'known indications',
        False: 'non-indications'
    })
    return scores


def rename_and_order_func(scores, nice_function_names):
    scores.func = scores.func.cat.rename_categories(nice_function_names)
    scores.func = Categorical(scores.func, ordered=True, categories=[
        name
        for name in nice_function_names.values()
        if name in scores.func.unique()
    ])
    return scores


def categorize_values(scores):
    categorical_columns = ['pert_iname', 'group', 'pert_idose', 'func', 'cell_id']
    for categorical in categorical_columns:
        if categorical in scores.columns:
            scores[categorical] = Categorical(scores[categorical])
    gc.collect()
    return scores


def process_scores(scores, names_map):
    scores = categorize_values(scores)
    scores = normalize_scores(scores, rescale=True, by_cell=False)
    scores = rename_and_order_func(scores, names_map)
    scores = add_helper_columns(scores)
    scores = rank_by_func(scores)
    return scores
