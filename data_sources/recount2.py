from pandas import read_table, concat

from config import DATA_DIR
from data_sources.data_source import DataSource
from data_sources.sra import SRAExpressionLayer
from layers import ExpressionLayer


class Recount2(DataSource):

    def __init__(self, meta_sra, *args, **kwargs):
        # if not meta_sra:
        super().__init__(*args, **kwargs)
        self.meta_sra = meta_sra

    def expression(self, studies):
        expressions = {}
        for study in studies:
            expressions[study] = read_table(
                DATA_DIR + f'/recount2/{study}.tsv.gz'
            ).set_index('gene_id')
        df = concat(
            expressions.values(),
            axis='columns'
        )
        layer = SRAExpressionLayer(df, run_to_study={
            run: study
            for run in df.columns
        })
        return self.with_disease_info(layer)

    def with_disease_info(self, layer):
        runs = layer.columns
        diseases = {
            run: self.meta_sra.run_to_disease(run)
            for run in runs
        }
        layer.run_to_disease = diseases
        return layer

    supported_layers = {
        ExpressionLayer: expression
    }
