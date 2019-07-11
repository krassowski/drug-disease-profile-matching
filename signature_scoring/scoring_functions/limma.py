from random import randint

from pandas import DataFrame
from rpy2.rinterface import RRuntimeError
from rpy2.robjects import r, globalenv
from rpy2.robjects import StrVector, ListVector
from rpy2.robjects.packages import importr

from data_sources.molecular_signatures_db import MolecularSignaturesDatabase
from enhanced_multiprocessing.cache_manager import multiprocess_cache_manager

from ..models.with_controls import ExpressionWithControls
from . import scoring_function
from .gsea import combine_gsea_results


db = MolecularSignaturesDatabase()
LIMMA_CACHE = None
multiprocess_cache_manager.add_cache(globals(), 'LIMMA_CACHE', 'dict')


def roast(expression: ExpressionWithControls, gene_sets: str, use_cache: bool):
    if use_cache:
        key = (expression.hashable, gene_sets)
        if key in LIMMA_CACHE:
            return LIMMA_CACHE[key]

    joined = DataFrame(expression.joined)
    joined.index = joined.index.astype(str)

    nulls = joined.isnull().any(axis=0).reset_index(drop=True)

    if nulls.any():
        print(f'Following columns contain nulls and will be skipped: {list(joined.columns[nulls])}')
        joined = joined[joined.columns[~nulls]]
        classes = expression.classes.loc[~nulls.reset_index(drop=True)]
    else:
        classes = expression.classes

    globalenv['expression'] = joined
    globalenv['expression_classes'] = classes

    # TODO: paired samples for cell lines?
    result = r(f"""
    expression_set = ExpressionSet(assayData=data.matrix(expression))

    design = cbind(intercept=1, controlVsCondition=expression_classes != 'normal')
    result = mroast(expression_set, {gene_sets}, design, geneid=dimnames(expression)[[1]])
    rows = rownames(result)
    result
    """)
    result.index = r['rows']

    r('rm(expression_set, design, expression, result, rows)')

    result = result.rename({'FDR': 'fdr_q-val'}, axis=1)
    result['nes'] = result[['PropUp', 'PropDown']].max(axis=1) * result.Direction.map({'Up': 1, 'Down': -1})
    if use_cache:
        LIMMA_CACHE[key] = result
    return result


def create_roast_scorer(
    gene_sets='c2.cp.kegg', id_type='entrez', grouping='by_substance',
    q_value_cutoff=0.1, na_action='fill_0', cache=True, cache_signatures=False
):
    """Only cache signatures when doing permutations, otherwise it will only slow it down"""

    importr('limma')
    importr('Biobase')

    gene_sets_r = ListVector({
        gene_set.name: StrVector(list(gene_set.genes))
        for gene_set in db.load(gene_sets=gene_sets, id_type=id_type).gene_sets
    })

    def set_gene_set_collection():
        globalenv[gene_sets] = gene_sets_r

    def roast_score(disease: ExpressionWithControls, compound: ExpressionWithControls):

        if len(compound.cases.columns) < 2 or len(compound.controls.columns) < 2:
            print(f'Skipping {compound} not enough degrees of freedom (no way to compute in-group variance)')
            return None

        if cache:
            multiprocess_cache_manager.respawn_cache_if_needed()

        try:
            disease_gene_sets = roast(disease, gene_sets=gene_sets, use_cache=cache)
            disease_gene_sets.drop(disease_gene_sets[disease_gene_sets['fdr_q-val'] > q_value_cutoff].index, inplace=True)

            signature_gene_sets = roast(compound, gene_sets=gene_sets, use_cache=cache and cache_signatures)

            joined = combine_gsea_results(disease_gene_sets, signature_gene_sets, na_action)

            if randint(0, 100) == 1:
                r('gc()')

            return joined.score.mean()
        except RRuntimeError as e:
            print(e)
            return None

    return scoring_function(
        roast_score, input=ExpressionWithControls, grouping=grouping,
        before_batch=set_gene_set_collection
    )
