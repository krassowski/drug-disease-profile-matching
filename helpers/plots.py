from copy import copy

from matplotlib import pyplot as plt
from matplotlib import patches
from pandas.core.dtypes.dtypes import CategoricalDtype
from plotnine.ggplot import ggplot


def y_line(p, score_mean, color='#bf4045', width=1.5, label='', n=1):
    x, y = p.get_lines()[n].get_data()
    y_by_dist_from_mean = {}
    for xi, yi in zip(x, y):
        dist = abs(xi - score_mean)
        y_by_dist_from_mean[dist] = yi
    closest = min(y_by_dist_from_mean.keys())

    plt.vlines(score_mean, 0, y_by_dist_from_mean[closest], linewidth=width, colors=color, label=label)


def pinecone_plot(p: ggplot, observed_column: str, maximized_rows=None, hatch='-------'):
    """Apply hatch to a part of given plot p, with the threshold determined by value from observed_column.
    maximized_rows: whether values over or below threshold should be hatched.

    In order to preserve correct ordering of rows/columns and categories,
    corresponding columns in data DataFrame have to be categorical.
    """

    maximized_rows = maximized_rows or []
    figure = p.draw()

    # suppress the figure drawn above from showing up
    plt.close()

    axes = figure.get_axes()

    columns_id = p.facet.cols[0]
    rows_id = p.facet.rows[0]

    data = p.data

    columns = data[columns_id]
    rows = data[rows_id]

    assert isinstance(columns.dtype, CategoricalDtype) and isinstance(rows.dtype, CategoricalDtype)

    columns = columns.unique()
    rows = rows.unique()

    categories = p.mapping['x']
    values = p.mapping['y']

    assert isinstance(data[categories].dtype, CategoricalDtype)

    axes_mapping = []
    for row in rows:
        for column in columns:
            axes_mapping.append((row, column))

    for i, a in enumerate(axes):
        row, column = axes_mapping[i]
        axis_data = data[data[columns_id] == column]
        axis_data = axis_data[axis_data[rows_id] == row]
        observed = axis_data[[observed_column, categories]].drop_duplicates().set_index(categories)
        collections = list(a.collections)

        for j, category in enumerate(observed.index):
            threshold = observed.loc[category][observed_column]
            collection = collections[j]
            category_values = axis_data[axis_data[categories] == category][values]
            y_max = category_values.max()
            y_min = category_values.min()

            maximized = row in maximized_rows
            start = threshold if maximized else y_min
            end = y_max - threshold if maximized else threshold - y_min

            hide_less_extreme = patches.Rectangle((0.5 + j, start), 1, end, transform=a.transData)

            new_c = a.add_collection(copy(collection))
            new_c.set_hatch(hatch)
            new_c.set_clip_path(hide_less_extreme)

    return figure
