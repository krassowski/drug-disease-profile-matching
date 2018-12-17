from enum import Enum
from types import SimpleNamespace, FunctionType
from typing import Dict, Tuple, List

from pandas import concat, DataFrame

from data_sources.drug_central import substances_for
from data_sources.drug_connectivity_map import dcm
from data_sources.census import CancerCensus
from utilities_namespace import show_table
from helpers.cache import cache_decorator
import metrics


class Mode(Enum):
    regression = 1
    classification = 2


census = CancerCensus().census


census_genes = [str(e) for e in census.entrez_geneid]
encoded_census_genes = [e.encode() for e in census_genes]


class AssociationsPreparation:

    def __init__(self, disease_code_to_disease_terms: Dict[str, List[str]], get_disease_expression: FunctionType):
        self.disease_code_to_disease_terms = disease_code_to_disease_terms
        controls_ids = dcm.all_controls(consensus=True).sig_id
        self.all_controls = dcm.from_ids(controls_ids, limit_to_one=True).reindex(encoded_census_genes)
        self.get_disease_expression = get_disease_expression

    # TODO: might need re-writing after change of get_differential (no longer ignores self by default; see ignore_first)
    @cache_decorator
    def get_differential(self, cancer_type, metric, limit_to_census, only_paired):

        expression_data = self.get_disease_expression(cancer_type, index='entrez_gene_id')
        differential = expression_data.differential(
            metric=metric,
            limit_to=census_genes if limit_to_census else None,
            only_paired=only_paired
        )
        return differential

    def find_substances(self, disease_name: str, contra: bool=False):
        search_terms = self.disease_code_to_disease_terms[disease_name]
        return [
            substance
            for term in search_terms
            for substance in substances_for(term, contra)
        ]

    def get_data_for(
        self,
        disease_name, limit_to_census=True, show_tables=False, ranks=False,
        mode=Mode.regression, cell_line_controls=True, disease_controls=True,
        fold_changes=True, only_paired=True
    ) -> Tuple[DataFrame, List]:

        assert not (fold_changes and ranks)
        assert mode in {mode.regression, mode.classification}

        metric = metrics.fold_change if fold_changes else metrics.signal_to_noise

        differential = self.get_differential(disease_name, metric, limit_to_census, only_paired)

        if differential is None:
            print(f'Not enough samples for {disease_name}, skipping')
            return None, None

        indicated_substances = self.find_substances(disease_name)
        contraindicated_substances = self.find_substances(disease_name, contra=True)

        indications = dcm.from_perturbations(indicated_substances, limit_to_one=True).reindex(encoded_census_genes)
        contraindications = dcm.from_perturbations(contraindicated_substances, limit_to_one=True).reindex(encoded_census_genes)

        controls = self.all_controls

        def fold_change(perturbation):
            control = dcm.get_control(perturbation.name, limit_to_genes=encoded_census_genes)
            assert (control.index == perturbation.index).all()
            return perturbation.divide(control)

        if fold_changes:
            indications = indications.apply(fold_change, axis='rows')
            contraindications = contraindications.apply(fold_change, axis='rows')
            controls = controls.apply(fold_change, axis='rows')

        common_genes = set(encoded_census_genes)
        common_genes -= set(indications[indications.isnull().any(axis=1)].index)
        common_genes -= set(contraindications[contraindications.isnull().any(axis=1)].index)
        common_genes -= set(differential[differential.isnull()].index)

        indications = indications.loc[common_genes]
        contraindications = contraindications.loc[common_genes]
        differential = differential.loc[common_genes]
        controls = controls.loc[common_genes]

        if ranks:
            indications = indications.rank(axis='rows')
            contraindications = contraindications.rank(axis='rows')
            differential = differential.rank()

        # drug inducing expression same as the disease (so basically causing a disease) is no cure for the disease
        disease_controls = DataFrame(differential, columns=[disease_name]) if disease_controls else DataFrame()

        differential.index = differential.index.map(lambda gene: b'disease_' + gene)

        input_data = []

        if not cell_line_controls:
            controls = controls.drop(controls.columns)

        categories = [indications, contraindications, controls, disease_controls]

        if show_tables:
            for category in categories:
                show_table(category, n_rows=3)

        for category in categories:
            signatures = category.columns

            if len(signatures):
                dummy = concat([differential] * len(signatures), axis=1)
                dummy.columns = signatures
                input_for_category = concat([category, dummy])
                input_data.append(input_for_category)

        if not input_data:
            return None, None

        x = concat(input_data, axis=1)

        labels_by_mode = {
            Mode.classification: SimpleNamespace(
                indication='indication',
                contraindication='contraindication',
                control='control'
            ),
            Mode.regression: SimpleNamespace(
                indication=1,
                contraindication=-1,
                control=0
            )
        }

        labels = labels_by_mode[mode]

        y = (
            [labels.indication] * len(indications.columns) +
            [labels.contraindication] * len(contraindications.columns) +
            [labels.control] * len(controls.columns) +
            [labels.control] * len(disease_controls.columns)
        )

        return x, y
