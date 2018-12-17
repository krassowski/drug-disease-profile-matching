from pandas import crosstab

from data_frames import AugmentedDataFrame
from helpers.cache import cached_property
from layers import Layer, CodingMutationLayer, MutationLayer


class LayerDataSubset:

    def __init__(self, filter, type):
        self.filter = filter
        self.type = type

    def value(self, layer_data):
        layer_data_type = layer_data.__class__

        return layer_data_type(
            self.filter(layer_data),
            layer_type=self.type
        )


class LayerData(AugmentedDataFrame):#, meta=ABCMeta):

    def as_layer(self, layer_type=None) -> Layer:
        pass

    @cached_property
    def layer(self) -> Layer:
        return self.as_layer()


class LayerDataWithSubsets(LayerData):
    # why sublcass, not a mixin? -> due to __getattr__ behaviour in mutli-inheritance

    __subsets__ = {}

    def __getattr__(self, key):
        try:
            subset = self.__subsets__[key]
            return subset.value(self)
        except KeyError:
            return super().__getattr__(key)


Subset = LayerDataSubset


class MutationAnnotations(LayerDataWithSubsets):
    """Subclass of DataFrame holding mutations data,

    which come from file in Mutation Annotation Format (MAF),
    adding layers generation support and utility functions.
    """

    __subsets__ = {
        'coding': Subset(
            filter=lambda mutations: (
                mutations[mutations.variant_classification.isin(
                    {'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Slient'}
                )]
            ),
            type=CodingMutationLayer
        ),
        'missense': Subset(
            filter=LayerData.having(variant_classification='Missense_Mutation'),
            type=CodingMutationLayer
        ),
        'nonsense': Subset(
            filter=LayerData.having(variant_classification='Nonsense_Mutation'),
            type=CodingMutationLayer
        ),
        'all_muts': Subset(
            filter=lambda x: x,
            type=MutationLayer
        )
    }

    def __init__(self, *args, layer_type=MutationLayer, **kwargs):
        super().__init__(*args, **kwargs)
        self.columns = [column.lower() for column in self.columns]
        self.layer_type = layer_type

    def propagate_counts(self):
        raise NotImplementedError

    def as_layer(self, layer_type=None):
        layer_type = layer_type if layer_type is not None else self.layer_type
        raw_layer = crosstab(self.participant, self.hugo_symbol)
        return layer_type(raw_layer)
