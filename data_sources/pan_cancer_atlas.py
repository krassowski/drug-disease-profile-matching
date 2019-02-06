from pandas import read_excel
from rpy2.robjects.packages import importr
from rpy2.robjects.pandas2ri import ri2py

from config import DATA_DIR
from data_frames import AugmentedDataFrame
from data_sources.data_source import DataSource
from data_sources.tcga import TCGA


class PanCancerAtlas(DataSource):

    url = 'https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/subtypes.html'

    def __init__(self):
        self.tcga_biolinks = importr('TCGAbiolinks')
        self.tcga = TCGA()

    def subtypes_curated(self):
        print(f'Please refer to {self.url} for the source of data for specific cohorts')
        pan_cancer_subtypes = self.tcga_biolinks.PanCancerAtlas_subtypes()
        pan_cancer_subtypes = ri2py(pan_cancer_subtypes)

        self.tcga.add_participant_column(pan_cancer_subtypes, 'pan.samplesID')

        return AugmentedDataFrame(pan_cancer_subtypes)

    def subtypes_integrative(self, tumor):
        pan_cancer_subtypes = self.tcga_biolinks.TCGAquery_subtype(tumor=tumor)
        pan_cancer_subtypes = ri2py(pan_cancer_subtypes)

        self.tcga.add_participant_column(pan_cancer_subtypes, 'patient')

        return pan_cancer_subtypes

    def icluster(self):
        df = read_excel(DATA_DIR + '/stratification/mmc6.xlsx', skiprows=1)

        # no information about samples - participants only
        assert (df['Sample ID'].str.len() == 12).all()

        df_with_participants = self.tcga.add_participant_column(df, 'Sample ID')

        # no duplicate participants
        assert not df_with_participants.participant.duplicated().any()

        return df_with_participants
