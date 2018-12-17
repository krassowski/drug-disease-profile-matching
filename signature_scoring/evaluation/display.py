from pandas import DataFrame, MultiIndex

from utilities_namespace import show_table
from helpers.gui.table import bordered_table
from helpers.gui.namespace import NeatNamespace
from . import metrics_manager


by_objective = metrics_manager.metrics_by_objective(convert_to=NeatNamespace)

minimized_metrics = by_objective.minimize.__dict__
maximized_metrics = by_objective.maximize.__dict__


def choose_columns(df: DataFrame, metrics_by_category: dict, categories: set):
    chosen_metrics = [
        f'{category}:{metric}'
        for category in categories
        for metric in metrics_by_category.get(category, [])
    ]
    return df.columns.intersection(chosen_metrics)


def highlight_best(column, color='#D2FAD2'):
    category, name = column.name
    return [
        f'background-color: {color}' if is_best else ''
        for is_best in metrics_manager.best_scores(category, name, scores=column)
    ]


def rank_and_sort(df, ranked_categories):
    maximized_columns = choose_columns(df, maximized_metrics, ranked_categories)
    minimized_columns = choose_columns(df, minimized_metrics, ranked_categories)

    df['overall:Rank score'] = (
        df[maximized_columns].rank().sum(axis=1)
        -
        df[minimized_columns].rank().sum(axis=1)
    )
    df = df.dropna(axis='columns')

    sort_columns = ['overall:Rank score', *maximized_columns, *minimized_columns]
    sort_ascending = [False] + [False] * len(maximized_columns) + [True] * len(minimized_columns)

    ranked_table = df.sort_values(sort_columns, ascending=sort_ascending)
    return ranked_table


def summarize(
    df,
    ranked_categories={'overall', 'indications', 'contraindications', 'controls'},
    sort='ranks',
    return_table=False
):
    ranked_table = rank_and_sort(df, ranked_categories)

    tuples = [
        column.split(':') if ':' in column else ('meta', column)
        for column in ranked_table.columns
    ]
    ranked_table.columns = MultiIndex.from_tuples(tuples, names=['Category', 'Metric'])

    display_columns = {
        column
        for column in ranked_table.columns
        if column[0] in ranked_categories or column[0] == 'meta' or column[1] == 'Rank score'
    }
    # I don't show the raw scores as this would inflate the table
    display_columns -= {('meta', 'Scores')}

    display_ready_table = ranked_table[list(display_columns)]
    display_ready_table = display_ready_table.sort_index(axis='columns')

    if sort != 'ranks':
        display_ready_table = display_ready_table.sort_values(sort, axis='rows', ascending=False)

    # this bit is used for fold-change benchmark results (see Benefits of additional fold-change...)
    if display_ready_table.index.duplicated().any():
        display_ready_table = display_ready_table.reset_index().set_index([('meta', 'Fold'), ('Func',)])

    show_table(
        display_ready_table
        .style
        .apply(highlight_best, axis=0)
        .set_table_styles(bordered_table(hide_headers=[3]))
        .background_gradient(cmap='RdYlGn', subset=[('overall', 'Rank score')], low=0.6, high=.8)
    )

    if return_table:
        return display_ready_table
