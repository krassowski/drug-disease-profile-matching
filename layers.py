from typing import Mapping, List

from pandas import DataFrame, Series, concat, merge
from sklearn.preprocessing import minmax_scale

from bio.protein_sequences import ProteinSequences
from data_frames import AugmentedDataFrame


class Layers:
    
    def __init__(self, data: Mapping[str, DataFrame]):
        self.layers = data

    def to_df(self) -> DataFrame:
        return concat([v for v in self.layers.values()], axis='columns')

    def rescale(self, scaler=minmax_scale) -> 'Layers':
        scaled = {}

        for key, value in self.layers.items():
            # scaled[key] = pd.DataFrame(value.apply(zscore, axis=0))
            scaled[key] = DataFrame(scaler(value))

        return Layers(scaled)

    @property
    def layers_with_negatives(self) -> List[str]:
        return [
            key
            for key, value in self.layers.items()
            if value.min().min() < 0
        ]

    def scale_proportionally_to_certainty(self, scaler=None):
        scaled = {}

        #from random import uniform
        #from numpy import nanmean
        from statistics import mean
        
        if scaler:
            self = self.rescale(scaler=scaler)

        from .nmf import MultiViewNMF
        # TODO: I use only single view features here, maybe another clas for that?
        nmf = MultiViewNMF()
        
        flip = True if self.layers_with_negatives else False
            
        for key, value in self.layers.items():

            # one could do supervised learning here, optimizing a scaling factor for each layer
            # or it could be related to the internal performance of clustering of individual layers
            # (if a layer does not know how to cluster into three parts, it is not very useful in this case)

            #if rescaled_layers.layers_with_negatives:

            n = nmf.flip_negatives(value) if flip else nmf.set_matrix(value)
            
            W, H = n.decompose()
            #p = W.predictions_certainty()
            p = W.predictions_certainty_combined()
            #p = W.first_prediction_certainty()
            #print(p)
            w = mean(p.fillna(0))
            print(key, w)


            # ww = value
            #while 1 - w > 0.04:
            #    ww = pd.DataFrame(ww) / w
            #    W, H = nmf.flip_negatives(ww).decompose()
            #    p = W.predictions_certainty()
            #    #p = W.predictions_certainty_combined()
            #    #p = W.first_prediction_certainty()
            #    w = mean(p)
            #    print(w)

            # here is an important observation: the more certain the layer,
            # the worse is the clustering. Why? Possibly there is a confounding
            # factor of the data types provided for each of the layers:
            # mutations are very sparse and single occurrences can lead to both:
            # incorrect clustering classification AND high certainty (if it is the only mutation)

            # dividing values by certainty leads to...
            # multiplicating values be certainty leads to..

            # maybe is is a wrong appraoch? maybe it has to to applied only for weighted,
            # combination of data from single layer clustering???? 

            # anyway, I could look into strategies for imputing/downgrading spurious mutaitons
            scaled[key] = value * w   # uniform(0, 1)

        return Layers(scaled)
    

class Layer(AugmentedDataFrame):
    """
    rows = samples, patients
    columns = genes, features
    """
    
    def limit_to_samples(self, samples_set):
        return self[self.index.isin(samples_set)]
    
    def rearange_to_match_with(self, other_layer, how='left', match='rows'):
        assert match == 'rows'
        index_only = other_layer[[]]
        rearranged = merge(
            index_only, self,
            left_on=index_only.index, right_on=self.index,
            how='left'
        ).set_index('key_0')
        return self._constructor(rearranged)
    
    def preprocess(self, rearrange_against=None):
        new_layer = (
            self.rearange_to_match_with(rearrange_against)
            if rearrange_against is not None else
            self
        )
        return (new_layer
            .drop_useless_features()
            .clean_index_names()
            .fillna(0)
        )
    
    def clean_index_names(self):
        self.columns.name = ''
        self.index.name = ''
        return self
    
    def drop_useless_features(self):
        without_na = self.dropna(axis='columns', how='all')
        #print()
        #print(without_na)
        unique_count = without_na.apply(Series.nunique)
        #print(unique_count)
        without_variance = unique_count[unique_count == 1].index
        #print(without_variance)
        return without_na.drop(without_variance, axis='columns')


class MutationLayer(Layer):
    
    __default_normalizer__ = None
    
    #@cached_property
    @property
    def normalized(self):
        return self.__default_normalizer__()

class ExpressionLayer(Layer):
    
    pass
    
# TODO pony database?

protein_sequences = ProteinSequences('data/uniprot/uniprot_sprot.fasta')

from numpy import nan

# if not hierarchy I can always refactor to mixins for more flexibility
class CodingMutationLayer(MutationLayer):

    def normalize_by_protein_length(self, recover=True):
        
        from statistics import StatisticsError
        failed = 0

        def normalize(column):
            try:
                protein = column.name
                protein_length = protein_sequences.average_length(protein)
                return column / protein_length
            except (StatisticsError, KeyError):
                nonlocal failed
                failed += 1
                if recover:
                    return column / protein_sequences.average_length()
                else:
                    return column * nan

        if failed:
            from warnings import warn
            warn(f'For {failed} proteins the exact lengths were not available; average protein length was used instead')

        return self.apply(normalize, axis='rows')
    
    __default_normalizer__ = normalize_by_protein_length

    
class NonCodingMutationLayer(MutationLayer):

    def normalize_by_gene_length(self):
        pass
