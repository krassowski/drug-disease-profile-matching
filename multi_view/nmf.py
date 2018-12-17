from pandas import DataFrame, concat
from copy import deepcopy
from sklearn.decomposition import NMF
from sklearn.metrics.cluster import adjusted_rand_score
from typing import Iterable


class DecomposedMatrix:
    
    def __init__(self, array, expected_clustering: DataFrame):
        self.array = array
        self.expected_clustering = expected_clustering
    
    def predict_clusters(self) -> DataFrame:
        clusters = DataFrame(self.array).idxmax(axis=1)
        clusters = clusters.to_frame()
        clusters.columns = ['cluster']
        return clusters
    
    def first_prediction_certainty(self):
        # [(cluster) = (weight)] => result 
        # [A = 1.0,  B = 0.5,  C = 0.5] => 1 - 0.5 = 0.5
        # [A = 1.0,  B = 0.5,  C = 0.0] => 1 - 0.5 = 0.5    # this one should get a higher score than 1,0.5,0.5
        # [A = 1.0,  B = 1.0,  C = 0.0] => 1 - 1.0 = 0
        # [A = 1.0,  B = 0.0,  C = 0.0] => 1 - 0.0 = 1
        m = DataFrame(self.array).max(axis=1)
        sm = DataFrame(self.array).apply(lambda row: sorted(row)[-2], axis=1)
        sum_ = DataFrame(self.array).apply(sum, axis=1)
        return m / sum_ - sm / sum_

    def predictions_certainty(self):
        # [(cluster) = (weight)] => result 
        # [A = 1.0,  B = 0.5,  C = 0.5] => 1 / (1 + 1) = 1 / 2 = 0.5
        # [A = 1.0,  B = 0.5,  C = 0.0] => 1 / (1 + 0.5) = 1 / 1.5 = 0.(6)
        # [A = 1.0,  B = 1.0,  C = 0.0] => 1 / (1 + 1) = 1 / 2 = 0.5   # this one should get a lower score than 1,0.5,0.5
        # [A = 1.0,  B = 0.0,  C = 0.0] => 1 / (1 + 0) = 1 / 1 = 1
        m = DataFrame(self.array).max(axis=1)
        sm = DataFrame(self.array).apply(sum, axis=1)
        return 1 / (1 + (sm - m)/m)

    def predictions_certainty_combined(self):
        return self.first_prediction_certainty() + self.predictions_certainty()
    
    def transform(self, function):
        new = deepcopy(self)
        new.array = function(new.array)
        return new
    
    def score(self, function=adjusted_rand_score):
        clusters = self.predict_clusters()
        return function(clusters.cluster, self.expected_clustering.cluster)


class MultiViewNMF:
    
    def __init__(self, matrix=None, expected_clustering: list=None, **kwargs):
        self.nmf = NMF(**kwargs)
        self.matrix = matrix
        self.expected_clustering = (
            self._clustering_frame(expected_clustering)
            if expected_clustering
            else None
        )

    def trim_negatives(self, matrix=None):
        matrix = matrix if matrix is not None else self.matrix
        assert matrix is not None
        trimmed = matrix.applymap(lambda x: max(x, 0))
        new = deepcopy(self)
        new.matrix = trimmed
        return new
        
    def _matrix_or_layers_or_own_matrix(self, matrix, layers):
        assert not (matrix is not None and layers is not None)
        if layers is not None:
            from layers import Layers
            if isinstance(layers, Layers):
                matrix = layers.to_df()
            else:
                # just dict
                matrix = concat(layers.values(), axis=1)
        matrix = matrix if matrix is not None else self.matrix
        assert matrix is not None
        return matrix

    def flip_negatives(self, matrix=None, layers=None):
        matrix = self._matrix_or_layers_or_own_matrix(matrix, layers)
        trimmed = matrix.applymap(lambda x: max(x, 0))
        reflected = matrix.applymap(lambda x: -min(x, 0))
        matrix = concat([trimmed, reflected], axis=1)
        new = deepcopy(self)
        new.matrix = matrix
        return new
    
    def set_matrix(self, matrix=None, layers=None):
        matrix = self._matrix_or_layers_or_own_matrix(matrix, layers)
        new = deepcopy(self)
        new.matrix = matrix
        return new
    
    def _clustering_frame(self, clustering: Iterable):
        return DataFrame(clustering, columns=['cluster'])

    def decompose(self, matrix=None, layers=None, exclude=None):
        matrix = self._matrix_or_layers_or_own_matrix(matrix, layers)
        expected_clustering = self.expected_clustering
        if exclude:
            rows_to_keep = ~matrix.index.isin(exclude)
            # matrix = copy(matrix)
            matrix = matrix[rows_to_keep]
            expected_clustering = expected_clustering[rows_to_keep]
        W = self.nmf.fit_transform(matrix)
        H = self.nmf.components_
        
        return DecomposedMatrix(W, expected_clustering), DecomposedMatrix(W, expected_clustering)
    
    def __repr__(self):
        return f'<MultiViewNMF of a matrix {self.matrix.shape}>'


class NIMFA(MultiViewNMF):
    
    #nsmf = factorization.nmf.Nmf(data.values, rank=3, max_iter=400, n_run=1)
    #nsmf.factorize()
    #classes = nsmf.predict(what='features')
    #return adjusted_rand_score(classes.tolist()[0], subtypes.subtype)
    pass