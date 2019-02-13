from copy import copy

from matplotlib import pyplot as plt, patches
from pandas.core.dtypes.dtypes import CategoricalDtype


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

    data = p.data

    try:
        columns_id = p.facet.cols[0]
        columns = data[columns_id]
        assert isinstance(columns.dtype, CategoricalDtype)
        columns = columns.unique()
    except IndexError:
        columns = [None]

    try:
        rows_id = p.facet.rows[0]
        rows = data[rows_id]
        assert isinstance(rows.dtype, CategoricalDtype)
        rows = rows.unique()
    except IndexError:
        rows = [None]

    categories = p.mapping['x']
    values = p.mapping['y']

    assert isinstance(data[categories].dtype, CategoricalDtype)

    axes_mapping = []
    for row in rows:
        for column in columns:
            axes_mapping.append((row, column))

    for i, a in enumerate(axes):
        row, column = axes_mapping[i]
        axis_data = data
        if column is not None:
            axis_data = axis_data[axis_data[columns_id] == column]
        if row is not None:
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
