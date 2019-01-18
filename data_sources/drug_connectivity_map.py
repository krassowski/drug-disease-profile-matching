from functools import lru_cache
import warnings

from pandas import read_table, DataFrame, Series, concat
from tqdm import tqdm

from config import DATA_DIR
from data_frames import MyDataFrame
from data_sources.data_source import DataSource
from h5py import File
from h5py.h5py_warnings import H5pyDeprecationWarning

from models import ExpressionProfile
from helpers.cache import cached_property


warnings.simplefilter("ignore", H5pyDeprecationWarning)


class PerturbationProfile(ExpressionProfile):

    def classes(self):
        pass


class DrugConnectivityMap(DataSource):

    def __init__(self):
        dataset_path = DATA_DIR + '/lincs/GSE92742'

        CMAP_PATH = dataset_path + '/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx'
        self.cmap_file = File(CMAP_PATH, mode='r')
        self.cmap = self.cmap_file['0']
        self.cell_info = read_table(dataset_path + '/GSE92742_Broad_LINCS_cell_info.txt.gz', low_memory=False)
        self.sig_metrics = read_table(dataset_path + '/GSE92742_Broad_LINCS_sig_metrics.txt.gz')
        self.sig_index_vector = self.meta['COL']['id'].value
        self.sig_index = {
            sig_id: i
            for i, sig_id in enumerate(self.sig_index_vector)
        }
        self.sig_info = read_table(dataset_path + '/GSE92742_Broad_LINCS_sig_info.txt.gz', low_memory=False)

    def signatures_treated_with(self, substance: str, pert_id=False):
        return self.metadata_for_perturbation(substance, pert_id=pert_id).sig_id

    def metadata_for_perturbation(self, substance: str, pert_id=False):
        sig_info = self.sig_info
        return sig_info[sig_info['pert_id' if pert_id else 'pert_iname'] == substance]

    def cell_for_perturbations(self, substances):
        all_cells = []
        for substance in substances:
            cells = self.metadata_for_perturbation(substance).cell_id
            all_cells.extend(cells)
        return all_cells

    @property
    def meta(self):
        return self.cmap['META']

    @property
    def matrix(self):
        return self.cmap['DATA']['0']['matrix']

    def signature_index(self, signature_id):
        try:
            if type(signature_id) is str:
                signature_id = signature_id.encode('utf-8')
            return self.sig_index[signature_id]
        except IndexError:
            print(
                f'Failed to find signature: {signature_id},'
                f'likely due to concurrency issues, retrying'
            )
            return self.signature_index(signature_id)

    def profile_by_signature(self, signature_id):
        column_index = self.signature_index(signature_id)
        return self.matrix[column_index]

    def profiles_by_signatures(self, signature_ids):
        indices_map = {
            signature_id: self.signature_index(signature_id)
            for signature_id in signature_ids
        }
        indices = sorted(indices_map.values())
        return self.matrix[indices], sorted(indices_map, key=indices_map.get)

    @cached_property
    def entrez_gene_ids(self):
        return self.meta['ROW']['id'].value

    def from_ids(self, signature_ids: Series, filter=True, **kwargs):
        if filter:
            signature_ids = self.filter_signatures(signature_ids, **kwargs)
        if not signature_ids:
            return DataFrame(index=self.entrez_gene_ids, columns=signature_ids)
        data, ordered_ids = self.profiles_by_signatures(signature_ids)
        return SignaturesData(
            data.T,
            index=self.entrez_gene_ids,
            columns=ordered_ids
        )

    def from_id(self, signature_id):
        index = self.signature_index(signature_id)
        data = self.matrix[index]
        return Series(data.T, index=self.entrez_gene_ids, name=signature_id)

    def iterate_signatures(self, signature_ids):
        for signature_id in signature_ids:
            index = self.signature_index(signature_id)
            data = self.matrix[index]
            yield SignaturesData(data.T, index=self.entrez_gene_ids, columns=[signature_id])

    def ids_of_exemplars(self, cell_id=None, take_first_per_pert=False):
        metrics = self.sig_metrics
        selected = metrics[metrics.is_exemplar.map(bool)]
        if cell_id:
            selected = selected.merge(self.sig_info)
            selected = selected[selected.cell_id == cell_id]
        if take_first_per_pert:
            selected = selected.drop_duplicates(subset=['pert_id'])
        return selected.sig_id

    def filter_signatures(self, signatures: Series, exemplar_only=True, limit_to_one=False, **kwargs):
        signatures = Series(signatures)
        if exemplar_only:
            metrics = self.sig_metrics
            signature_metrics = metrics[metrics.sig_id.isin(signatures)]
            exemplar_signature_metrics = signature_metrics[signature_metrics.is_exemplar.map(bool)]
            if not exemplar_signature_metrics.empty:
                signatures = exemplar_signature_metrics['sig_id']

        for key, value in kwargs.items():
            if value is None:
                continue
            signatures_data = self.sig_info[self.sig_info.sig_id.isin(signatures)]
            signatures = signatures_data[signatures_data[key] == value]['sig_id']

        if limit_to_one:
            signatures_data = self.sig_info[self.sig_info.sig_id.isin(signatures)]
            signatures_data = signatures_data.sort_values(
                ['pert_idose', 'pert_time', 'cell_id'],
                ascending=[False, False, False]
            )
            signatures = list(signatures_data.sig_id)
            signatures = Series(signatures[:1])
        return signatures.tolist()

    def ids_for_perturbations(self, substances, synonyms_source=None, pert_id=False, **kwargs):
        chosen_signatures = []
        for substance in substances:
            names = [substance]
            if synonyms_source:
                synonyms = synonyms_source(substance)
                names.extend(synonyms)
            while names:
                name = names.pop(0)
                signatures = self.signatures_treated_with(name, pert_id=pert_id)
                signatures = self.filter_signatures(signatures, **kwargs)
                if signatures:
                    if name != substance:
                        print(f'Didn\'t find anything for {substance}, but got matches for synonym {name}.')
                    chosen_signatures.extend(signatures)
                    break
            #if len(signatures) == 0:
            #    logging.log(f'No signatures id for substance {substance}')
        return chosen_signatures

    def from_perturbations(self, substances, synonyms_source=None, exemplar_only=True, cell_id=None, limit_to_one=False, pert_id=False):
        chosen_signatures = self.ids_for_perturbations(
            substances, synonyms_source=synonyms_source,
            exemplar_only=exemplar_only, limit_to_one=limit_to_one, cell_id=cell_id,
            pert_id=pert_id
        )
        chosen_signatures = list(set(chosen_signatures))
        matched = set(self.identify_substances(chosen_signatures))
        substances = set(substances)
        diff = matched - substances
        if diff:
            if synonyms_source:
                print(f'Following substances were inferred indirectly: {diff}')
            elif pert_id:
                print(f'Following substances have names different than ids: {diff}')
            else:
                assert diff == set()

        no_matches = f'No matches for: {substances - matched}.' if substances - matched else ''

        print(
            f'Got perturbations for {len(matched)}/{len(substances)} or '
            f'{(len(matched)/len(substances) if len(substances) else 0)*100:.2f}% '
            f'of substances. {no_matches}'
        )
        return self.from_ids(chosen_signatures, filter=False)

    cache = {}

    def controls_untreated(self, consensus):
        subtype = '.cns' if consensus else ''
        return self.sig_info[self.sig_info.pert_type == f'ctl_untrt{subtype}']

    def controls_vehicle(self, consensus):
        subtype = '.cns' if consensus else ''
        return self.sig_info[self.sig_info.pert_type == f'ctl_vehicle{subtype}']

    def controls_vector(self, consensus):
        subtype = '.cns' if consensus else ''
        return self.sig_info[self.sig_info.pert_type == f'ctl_vector{subtype}']

    def all_controls(self, consensus):
        controls = [
            self.controls_untreated,
            self.controls_vector,
            self.controls_vehicle
        ]
        return concat([get_controls(consensus=consensus) for get_controls in controls])

    def get_controls(self, signature_id, limit_to_genes=None, exemplar_only=False):
        """
        Documentation from LINCS states that DMSO is the control for compound treatments,
        while empty vectors and other non-gene-coding inserts (e.g LacZ) are controls
        for genetic perturbagens.
        """
        sig_data = self.sig_info[self.sig_info.sig_id == signature_id].squeeze()
        cell_id = sig_data.cell_id
        pert_itime = sig_data.pert_itime

        if sig_data.pert_type == 'trt_cp':  # compound treatment
            pert_type = 'ctl_vehicle'
            pert_iname = 'DMSO'
        else:
            # TODO: how should I choose adequate vector?
            pert_type = 'ctl_vector'
            pert_iname = None
            assert False

        signatures = self.signatures_treated_with(pert_iname)

        signatures = self.filter_signatures(
            signatures, pert_type=pert_type, pert_itime=pert_itime,
            cell_id=cell_id, exemplar_only=exemplar_only
        )

        controls = self.from_ids(signatures, filter=False)

        if limit_to_genes is not None:
            controls = controls.reindex(limit_to_genes)

        return controls

    def identify_substances(self, signature_ids):
        ids = set(signature_ids)
        return set(self.sig_info[self.sig_info.sig_id.isin(ids)].pert_iname)


from collections import UserDict
from statistics import mean


dcm = DrugConnectivityMap()


class AggregatedScores(DataFrame):
    pass


class Scores(UserDict):

    def __init__(self, *args, scores_for='sig_id', **kwargs):
        super().__init__(*args, **kwargs)
        self.df = DataFrame.from_dict(dict(self), orient='index', columns=['score'])
        self.merged = self.df.merge(
            dcm.sig_info,
            left_on=self.df.index, right_on=dcm.sig_info[scores_for]
        ).drop(['key_0', 'distil_id'], axis='columns')

    @classmethod
    def from_grouped_signatures(cls, data):
        per_single_signature = {}
        for signature_ids, score in data:
            for signature_id in signature_ids:
                if not score:
                    continue
                per_single_signature[signature_id] = score
        return cls(per_single_signature)

    def __add__(self, other):
        return Scores({**self, **other})

    @property
    def best_per_substance(self) -> AggregatedScores:
        return AggregatedScores(self.merged.groupby('pert_iname')['score'].max().sort_values(ascending=False))

    def limit_to_cell_line(self, cell_id):
        signature_ids = set(self.merged[self.merged.cell_id == cell_id].sig_id)
        new_scores = Scores(self.df.loc[signature_ids].score.to_dict())
        return new_scores

    def aggregate(self, by, func_name):
        if self.merged.empty:
            return AggregatedScores(columns=['score'])
        grouped = self.merged.groupby(by)['score']
        aggregated = getattr(grouped, func_name)()
        return AggregatedScores(aggregated.sort_values(ascending=False))

    @property
    def mean_per_substance_and_dose(self) -> AggregatedScores:
        return self.aggregate(['pert_iname', 'pert_idose'], 'mean')

    @property
    def mean_per_substance_dose_and_cell(self) -> AggregatedScores:
        return self.aggregate(['pert_iname', 'pert_idose', 'cell_id'], 'mean')

    @property
    def mean_per_substance(self) -> AggregatedScores:
        return self.aggregate('pert_iname', 'mean')

    @property
    def median_per_substance(self) -> AggregatedScores:
        return self.aggregate('pert_iname', 'median')

    @property
    def signal_to_noise(self):
        """
        Select only substances with more than one replicate.
        divide mean value by variation
        """
        grouped = self.merged.groupby('pert_iname')
        grouped = grouped.filter(lambda group: group['score'].count() > 1)
        grouped = grouped.groupby('pert_iname')
        if self.merged.empty:
            return Series()
        corrected = grouped['score'].mean() / grouped['score'].var()
        nans = corrected[corrected.isnull()]
        if not nans.empty:
            print('Nans detected, dropping', nans)
            corrected = corrected.dropna()
        return DataFrame(corrected.sort_values(ascending=False))

    @property
    def mean(self):
        return mean(list(self.values()))

    def _repr_html_(self):
        return self.merged.sort_values('score', ascending=False)._repr_html_()


class SignaturesData(MyDataFrame):

    def differential(self, dcm: DrugConnectivityMap, metric='difference_of_means'):
        tqdm.pandas()

        if metric == 'difference_of_means':
            def diff(column):
                controls = dcm.get_controls(column.name)
                mean_control = controls.mean(axis=1)
                return column - mean_control
        else:
            def diff(column):
                controls = dcm.get_controls(column.name)
                mean_control = controls.mean(axis=1)
                return (column - mean_control) / controls.std(axis=1)

        return self.progress_apply(diff, axis=0)

    def classes(self, class_type='pert_iname'):
        metadata = dcm.sig_info[dcm.sig_info.sig_id.isin(self.columns)]
        return metadata[class_type]

    @property
    def metadata(self):
        return dcm.sig_info[dcm.sig_info.sig_id.isin(self.columns)]

    def members_of_class(self, class_name, class_type='pert_iname', cell_line=None):
        metadata = self.metadata
        if cell_line:
            metadata = metadata[metadata['cell_id'] == cell_line]
        return metadata[metadata[class_type] == class_name].sig_id


dcm = DrugConnectivityMap()


@lru_cache()
def get_controls_for_signatures(ids, genes_to_keep=None):
    controls_by_signature = {}
    for signature_id in ids:
        controls = dcm.get_controls(signature_id, exemplar_only=True)
        if genes_to_keep is not None:
            rows_to_keep = controls.index.isin(genes_to_keep)
            controls = controls[rows_to_keep]
        control = controls.mean(axis=1)
        controls_by_signature[signature_id] = control
    return DataFrame(controls_by_signature)
