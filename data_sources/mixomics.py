from helpers.r import importr, data as package_data, data_frame_from_matrix

from .data_source import DataSource


class MixOmics(DataSource):
    """
    Installation:

        >>> from rpy2.robjects.packages import importr
        >>> utils = importr('utils')
        >>> utils.install_packages('mixOmics')
        >>> # installed 21 October, then update 6 November
    """

    def __init__(self):
        mixomics = importr('mixOmics')
        self.datasets = package_data(package=mixomics)

    def get_raw_dataset(self, name, subset):
        breast = self.datasets.fetch(name)
        data = breast[name]
        data_dict = dict(data.items())
        data = dict(data_dict[subset].items())
        return data
    
    def get_dataset(self, name='breast.TCGA', subset='data.train'):
    
        data = self.get_raw_dataset(name, subset)
        subtype = data.pop('subtype')
        subtype = list(subtype)
    
        data = {k: data_frame_from_matrix(v) for k, v in data.items()}
        
        return data, subtype
