from pandas import Series

from better_abc import abstract_property
from data_frames import AugmentedDataFrame


class ExpressionProfile(AugmentedDataFrame):

    @abstract_property
    def classes(self) -> Series:
        pass
