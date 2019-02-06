from functools import partial
from tempfile import NamedTemporaryFile
from typing import Union, Set
from warnings import warn

from pandas import concat

from data_frames import AugmentedDataFrame
from data_sources.molecular_signatures_db import MolecularSignaturesDatabase
from helpers.temp import create_tmp_dir
from methods.gsea import cudaGSEA
from methods.gsea.base import GSEA
from methods.gsea.exceptions import GSEAError, GSEANoResults
from methods.gsea.java import GSEADesktop
from multiprocess.cache_manager import multiprocess_cache_manager

from ..models import Profile
from ..models.with_controls import ExpressionWithControls, DummyExpressionsWithControls
from . import scoring_function, ScoringError

GSEA_CACHE = None
multiprocess_cache_manager.add_cache(globals(), 'GSEA_CACHE', 'dict')

tmp_dir = create_tmp_dir('gsea')


def combine_gsea_results(disease_gene_sets, signature_gene_sets, na_action='fill_0'):

    joined = disease_gene_sets.merge(
        signature_gene_sets, suffixes=['_disease', '_signature'],
        left_on=disease_gene_sets.index, right_on=signature_gene_sets.index,
        how=('left' if na_action == 'fill_0' else 'inner')
    )
    joined = joined.set_index('key_0')

    if na_action == 'fill_0':
        joined = joined.fillna(0)

    if (joined.nes_disease == 0).any():
        warn('ES in disease is 0, this will lead to inf scores!')

    joined['score'] = (
        1
        -
        (joined.nes_disease + joined.nes_signature) / joined.nes_disease
        *
        (1 - joined['fdr_q-val_disease']) * (1 - joined['fdr_q-val_disease'])
    )

    return joined


def cached_gsea_run(
    gsea_app,
    gsea, gene_sets, expression: ExpressionWithControls,
    class_name, warn_when_not_using_cache=False, delete=True
):

    profile_hash = hash(expression.hashable)
    key = (gene_sets, profile_hash, class_name, gsea_app.__class__.__name__)

    if key in GSEA_CACHE:
        results = GSEA_CACHE[key]
        if results is None:
            raise GSEANoResults()
    else:
        if warn_when_not_using_cache:
            warn('Warning: not using cache (that\'s fine if it\'s the first run)')
        results = gsea(
            expression,
            out_dir=f'{tmp_dir}/{class_name}',
            name=str(profile_hash).replace('-', 'm'),
            delete=delete
        )
        # if results is not None:
        GSEA_CACHE[key] = results

    return results[expression.case_name], results[expression.control_name]


DEFAULT_GSEA_APP = GSEADesktop()


def create_gsea_scorer(
    gsea_app: GSEA = DEFAULT_GSEA_APP, permutations=500, gene_sets='h.all',
    q_value_cutoff=0.1, na_action='fill_0', score='mean',
    normalization=True, metric='Diff_of_Classes',
    permutation_type='Gene_set', grouping=None,
    custom_multiprocessing=False, verbose=False,
    min_genes=15, max_genes=500, id_type='entrez',
    genes: Set[str] = None
):
    """
    na_action: fill_0 or drop
    score: mean, max, sum
        # mean = mean improvement for the condition (balancing pros and cons)
        # max = best effect on a molecular pathway, though might lead to deterioration of certain pathways of the disease.
        # sum =
    """

    molecular_signatures_db = MolecularSignaturesDatabase()

    with NamedTemporaryFile(delete=False, dir=tmp_dir, suffix='.gmt') as f:
        gene_sets_path = f.name
        matrix = molecular_signatures_db.load(gene_sets, id_type)
        before = len(matrix.gene_sets)

        if genes:
            matrix = matrix.subset(genes)

        matrix = matrix.trim(min_genes, max_genes)

        if verbose:
            print(f'Trimmed gene sets database from {before} to {len(matrix.gene_sets)}')
        matrix.to_gmt(gene_sets_path)

    if isinstance(gsea_app, cudaGSEA) and not genes and (min_genes or max_genes):
        warn(
            'Please supplement list of genes on the expression matrix '
            'to enable correct trimming of the gene sets for cudaGSEA'
        )

    input = Profile if permutation_type != 'phenotype' else ExpressionWithControls
    assert permutation_type in {'Gene_set', 'phenotype'}

    def to_dummy_expression(profile: Profile, class_name: str):
        case_expression = AugmentedDataFrame(
            concat([profile.top.up, profile.top.down])
        )
        return DummyExpressionsWithControls.from_differential(case_expression, case_name=class_name)

    def gsea_score(
        disease: Union[Profile, ExpressionWithControls],
        compound: Union[Profile, ExpressionWithControls],
        warn_about_cache=True
    ):
        # theoretically this could be replaced with the original data (e.g. tumour/normal for TCGA),
        # though this would be less consistent with how the signature profile is handled and require
        # re-calculation of the differential profile for each GSEA run.
        multiprocess_cache_manager.respawn_cache_if_needed()

        if isinstance(disease, Profile):
            disease_expression = to_dummy_expression(disease, 'disease')
            compound_expression = to_dummy_expression(compound, 'signature')
        else:
            disease_expression = disease
            compound_expression = compound

        gsea = partial(
            gsea_app.run,
            gene_sets=gene_sets_path, id_type=id_type,
            permutations=permutations, permutation_type=permutation_type,
            metric=metric,
            normalization='meandiv' if normalization else None,
            verbose=verbose,
            min_genes=min_genes, max_genes=max_genes
        )

        try:
            disease_gene_sets_up, disease_gene_sets_dn = cached_gsea_run(
                gsea_app,
                gsea, gene_sets, disease_expression, class_name='disease',
                warn_when_not_using_cache=warn_about_cache
                # delete=False might be beneficial for single runs (if these were to be restarted after the cache is gone)
                # but would also fill the disk with permutations quickly
            )
        except GSEAError as e:
            raise ScoringError(f'Diseases scoring failed, all the substances are doomed to fail; original error: {e}')

        try:
            signature_gene_sets_up, signature_gene_sets_dn = cached_gsea_run(
                gsea_app,
                gsea, gene_sets, compound_expression, class_name='signature'
            )
        except GSEAError:
            return None

        results = [disease_gene_sets_up, disease_gene_sets_dn]

        if q_value_cutoff:
            for result in results:
                result.drop(result[result['fdr_q-val'] > q_value_cutoff].index, inplace=True)

        disease_gene_sets = concat([disease_gene_sets_up, disease_gene_sets_dn])
        signature_gene_sets = concat([signature_gene_sets_up, signature_gene_sets_dn])

        joined = combine_gsea_results(disease_gene_sets, signature_gene_sets, na_action)

        return getattr(joined.score, score)()

    gsea_score.metadata = {
        'app': gsea_app.__class__.__name__,
        'permutations': permutations,
        'gene_sets': gene_sets,
        'score_calculation': score,
        'q_value_cutoff': q_value_cutoff,
        'na_action': na_action,
        'normalization': normalization,
        'permutation_type': permutation_type,
        'metric': metric
    }

    gsea_score.__name__ = '_'.join(map(str, gsea_score.metadata.values()))

    return scoring_function(
        gsea_score, input=input, grouping=grouping,
        custom_multiprocessing=custom_multiprocessing,
        before_batch=lambda: gsea_app.prepare_output(),
        supports_cache=True
    )

