from pandas import DataFrame, Series

from helpers.gui.namespace import HorizontalNamespace
from helpers.cache import power_hash


class AugmentedDataFrame(DataFrame):

    def normalize_columns(self):
        self.columns = [c.lower().replace(' ', '_') for c in self.columns]

    def having(self=None, **kwargs):

        class_ = __class__
        if self is not None:
            class_ = self.__class__

        def having_(new_layer):
            for key, value in kwargs.items():
                column = getattr(new_layer, key)
                new_layer = new_layer[column == value]
            return class_(new_layer)

        if self is not None:
            return having_(self)
        else:
            return having_

    @property
    def summary(self):
        data = {}
        for column in self.columns:
            column_summary = DataFrame(self[column].value_counts())
            if len(column_summary) == len(self):
                column_summary = f'{len(column_summary)} distinct records'
            data[column] = column_summary
        return HorizontalNamespace(data)

    def __hash__(self):
        return hash(power_hash(self))

    @property
    def _constructor(self):
        return self.__class__


class AugmentedSeries(Series):

    @property
    def _constructor(self):
        return self.__class__


def explode(df: DataFrame, column: str, keep_order=False):
    data = []
    columns = list(df.columns)
    to_explode_index = columns.index(column)
    columns.remove(column)
    columns = [column, *columns]
    for row in df.itertuples(index=False):
        base = [*row[:to_explode_index], *row[to_explode_index + 1:]]
        for entry in row[to_explode_index]:
            data.append([entry, *base])
    del df
    new_df = DataFrame(data, columns=columns)
    if keep_order:
        new_df = new_df[columns]
    return new_df


MyDataFrame = AugmentedDataFrame
