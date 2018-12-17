import sys

from pandas import DataFrame, Series

from data_sources.drug_connectivity_map import Scores, dcm
from helpers import first, WarningManager

from multiprocess import Pool
from multiprocess.cache_manager import multiprocess_cache_manager
from signature_scoring.profile import Signature, Profile

CACHE = None
multiprocess_cache_manager.add_cache(globals(), 'CACHE', 'dict')


def dummy_tqdm(a, *args, **kwargs):
    return a


progress_bar = dummy_tqdm


def switch_progress_bar(progress):
    global progress_bar
    if progress:
        from tqdm import tqdm
        progress_bar = tqdm
    else:
        progress_bar = dummy_tqdm


def select_genes(common_genes, gene_subset):
    if not gene_subset:
        return common_genes

    if type(first(common_genes)) is bytes and type(first(gene_subset)) is not bytes:
        gene_subset = {str(gene).encode() for gene in gene_subset}

    return common_genes & set(gene_subset)


class SignatureProcessor:

    signature_type = Signature

    def __init__(self, signatures: DataFrame, warning_manager=None, progress=False, processes=None):

        if signatures is None:
            print('No signatures given: using one exemplar signature per perturbagen from VCAP cell line')
            ids = dcm.ids_of_exemplars(cell_id='VCAP', take_first_per_pert=True).tolist()
            signatures = dcm.from_ids(ids[:1])
        else:
            ids = signatures.columns

        switch_progress_bar(progress)

        self.ids = ids
        self.progress = progress
        self.signatures = signatures
        self.processes = processes
        self.warning_manager = warning_manager or WarningManager()
        self.scale = False

    def get_signature(self, signature_id):
        if signature_id in self.signatures.columns:
            signature = self.signatures[signature_id]
        else:
            if signature_id in CACHE:
                signature = CACHE[signature_id]
            else:
                signature = dcm.from_id(signature_id)
                CACHE[signature_id] = signature
        return signature

    def score_signature(
        self, signature_id, disease_profile, rows_of_selected_genes, limit, scoring_func, gene_selection
    ):
        signature = self.get_signature(signature_id)
        signature = signature[rows_of_selected_genes]
        signature = self.transform_signature(signature, signature.index)
        signature = self.signature_type(signature)

        compound_profile = Profile(signature, limit, nlargest=gene_selection)

        score = scoring_func(disease_profile, compound_profile)

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

        common_genes = set(self.signatures.index) & set(disease_signature.index)
        common_genes = select_genes(common_genes, gene_subset)

        n = len(common_genes)

        self.warning_manager.warn_once(
            f'Retaining {n} genes: {n / len(self.signatures.index) * 100:.2f}% of signature genes and'
            f' {n / len(disease_signature.index) * 100:.2f}% of query genes'
        )

        return common_genes

    def warn_if_few_genes_selected(self, selected_genes, limit):
        if len(selected_genes) != 2 * limit and limit != sys.maxsize:
            self.warning_manager.warn_once(f'Selected only {len(selected_genes)} genes out of {2 * limit} allowed.')

    def score_signatures(
        self, scoring_func, disease_signature, limit=500, gene_subset=None,
        scale=False, gene_selection=Series.nlargest,
    ):
        # scaling will be performed in transform_signature
        self.scale = scale

        limit = limit or sys.maxsize

        common_genes = self.select_common_genes(disease_signature, gene_subset)
        disease_signature = disease_signature[disease_signature.index.isin(common_genes)]
        disease_signature = self.signature_type(disease_signature)

        disease_profile = Profile(disease_signature, limit=limit)

        selected_genes = disease_profile.top.genes
        self.warn_if_few_genes_selected(selected_genes, limit)

        rows_of_selected_genes = self.signatures.index.isin(selected_genes)

        shared_args = [disease_profile, rows_of_selected_genes, limit, scoring_func, gene_selection]

        # first signature is scored in one process,
        # so that the common cache is populated
        # without repetition of calculations
        scores = [self.score_signature(self.ids[0], *shared_args)]

        # and then iteratively apply scoring function to each next compound signature
        scores.extend(
            self.pool.imap(
                self.score_signature,
                self.ids[1:],
                shared_args=shared_args
            )
        )

        return Scores(scores)
