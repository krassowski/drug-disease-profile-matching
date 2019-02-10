import sys

from pandas import Series
from tqdm import tqdm

from data_sources.drug_connectivity_map import Scores, dcm
from helpers import first, WarningManager

from multiprocess import Pool
from multiprocess.cache_manager import multiprocess_cache_manager

from ..models import Signature, Profile, SignaturesGrouping
from ..scoring_functions import ScoringFunction


CACHE = None
multiprocess_cache_manager.add_cache(globals(), 'CACHE', 'dict')


def select_genes(common_genes, gene_subset):
    if not gene_subset:
        return common_genes

    if type(first(common_genes)) is bytes and type(first(gene_subset)) is not bytes:
        gene_subset = {str(gene).encode() for gene in gene_subset}

    return common_genes & set(gene_subset)


class SignatureProcessor:
    # TODO: this class in a somehow transitional state, between being a single-signature
    #  data wrangler and a data wrangler that supports grouped signatures (i.e. replicates/repetitions),
    #  therefore some of the naming might be a bit inconsistent

    signature_type = Signature
    scores_type = Scores

    def __init__(self, signatures: SignaturesGrouping, warning_manager=None, progress=False, processes=None):

        if signatures is None:
            print('No signatures given: using one exemplar signature per perturbagen from VCAP cell line')
            ids = dcm.ids_of_exemplars(cell_id='VCAP', take_first_per_pert=True).tolist()
            signatures = dcm.from_ids(ids[:1])
        else:
            ids = signatures.groups_keys()

        self.ids = ids
        self.progress = progress
        self.signature_groups: SignaturesGrouping = signatures
        self.processes = processes
        self.warning_manager = warning_manager or WarningManager()
        self.scale = False

    def get_signature_group(self, group_id):
        if group_id in self.signature_groups.groups_keys():
            signature = self.signature_groups[group_id]
        else:
            if group_id in CACHE:
                signature = CACHE[group_id]
            else:
                signature = dcm.from_id(group_id)
                CACHE[group_id] = signature
        return signature

    def score_signature_group(
        self, signature_id, disease_profile, rows_of_selected_genes, limit,
        scoring_func: ScoringFunction, gene_selection,
        warn_about_cache=True
    ):
        signature = self.get_signature_group(signature_id)
        signature = signature[rows_of_selected_genes]
        signature = self.transform_signature(signature, signature.index)

        if scoring_func.input == Profile:
            compound_profile = Profile(
                self.signature_type(signature),
                limit, nlargest=gene_selection
            )
        else:
            compound_profile = signature
            # TODO: apply limit to the compound_profile for the ExpressionsWithControls case.
            #  How? One idea: take the means/medians of gene values and choose n best genes.

        args = {}

        if scoring_func.custom_multiprocessing:
            args['cores'] = self.processes

        if scoring_func.supports_cache:
            args['warn_about_cache'] = warn_about_cache

        score = scoring_func(disease_profile, compound_profile, **args)

        del signature, compound_profile

        return signature_id, score

    @property
    def pool(self):
        return Pool(self.processes, progress_bar=self.progress)

    def transform_signature(self, signature, selected_genes):
        if self.scale:
            signature = signature / (signature.max() - signature.min())
        return signature

    def select_common_genes(self, disease_signature, gene_subset=None):

        common_genes = set(self.signature_groups.genes) & set(disease_signature.index)
        common_genes = select_genes(common_genes, gene_subset)

        n = len(common_genes)

        self.warning_manager.warn_once(
            f'Retaining {n} genes: {n / len(self.signature_groups.genes) * 100:.2f}% of signature genes and'
            f' {n / len(disease_signature.index) * 100:.2f}% of query genes'
        )

        return common_genes

    def warn_if_few_genes_selected(self, selected_genes, limit):
        if len(selected_genes) != 2 * limit and limit != sys.maxsize:
            self.warning_manager.warn_once(f'Selected only {len(selected_genes)} genes out of {2 * limit} allowed.')

    def single_process_map_with_shared(self, func, iterable, shared_args):
        if self.progress:
            iterable = tqdm(iterable)
        return [func(i, *shared_args) for i in iterable]

    def score_signatures(
        self, scoring_func, disease_signature, limit=500, gene_subset=None,
        scale=False, gene_selection=Series.nlargest, force_multiprocess_all=False
    ):
        # scaling will be performed in transform_signature
        self.scale = scale

        limit = limit or sys.maxsize

        common_genes = self.select_common_genes(disease_signature, gene_subset)
        disease_signature = disease_signature[disease_signature.index.isin(common_genes)]

        if scoring_func.input == Profile:
            disease_profile = Profile(
                self.signature_type(disease_signature),
                limit=limit, nlargest=gene_selection
            )
            selected_genes = disease_profile.top.genes
        else:
            disease_profile = disease_signature
            selected_genes = disease_profile.index

        self.warn_if_few_genes_selected(selected_genes, limit)

        rows_of_selected_genes = self.signature_groups.genes.isin(selected_genes)

        shared_args = [disease_profile, rows_of_selected_genes, limit, scoring_func, gene_selection]

        start = 0
        scores = []

        if scoring_func.multiprocessing_exclude_first and not force_multiprocess_all:
            # first signature is scored in one process,
            # so that the common cache is populated
            # without repetition of calculations
            scores = [self.score_signature_group(self.ids[0], *shared_args)]
            start = 1

        map_with_shared = (
            self.single_process_map_with_shared
            if scoring_func.custom_multiprocessing else
            self.pool.imap
        )

        # and then iteratively apply scoring function to each next compound signature
        scores.extend(
            map_with_shared(
                self.score_signature_group,
                self.ids[start:],
                shared_args=shared_args
            )
        )

        scores = [
            (signature_id, score)
            for signature_id, score in scores
            if score is not None
        ]

        if scoring_func.grouping:
            return Scores.from_grouped_signatures(scores)
        return Scores(scores)
