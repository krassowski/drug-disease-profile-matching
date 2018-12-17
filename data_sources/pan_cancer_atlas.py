from rpy2.robjects.packages import importr
from rpy2.robjects.pandas2ri import ri2py

from data_frames import AugmentedDataFrame
from data_sources.data_source import DataSource
from data_sources.tcga import TCGA


class PanCancerAtlas(DataSource):

    def __init__(self):
        self.tcga_biolinks = importr('TCGAbiolinks')
        self.tcga = TCGA()

    def subtypes(self):
        pan_cancer_subtypes = self.tcga_biolinks.PanCancerAtlas_subtypes()
        pan_cancer_subtypes = ri2py(pan_cancer_subtypes)

        self.tcga.add_participant_column(pan_cancer_subtypes, 'pan.samplesID')

        return AugmentedDataFrame(pan_cancer_subtypes)
