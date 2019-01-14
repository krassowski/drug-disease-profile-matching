from pathlib import Path
from typing import Union
from warnings import warn

from pandas import DataFrame, concat, Series
from rpy2.robjects import r, globalenv
from rpy2.robjects.packages import importr

from methods.gsea import MolecularSignaturesDatabase
from multiprocess.cache_manager import multiprocess_cache_manager

from ..models import ExpressionWithControls, Profile
from . import scoring_function
from .gsea import combine_gsea_results


db = MolecularSignaturesDatabase()
GSVA_CACHE = None
multiprocess_cache_manager.add_cache(globals(), 'GSVA_CACHE', 'dict')


process_specific_counter = 0


def gsva(expression: Union[ExpressionWithControls, Profile], gene_sets: str, method: str = 'gsva', single_sample=False, permutations=1000, mx_diff=True, _cache=True):
    """
    Excerpt from GSVA documentation:
        An important argument of the gsva() function is the flag mx.diff which is set to TRUE by default.

        Under this default setting, GSVA enrichment scores are calculated using Equation 5, and therefore, are
        more amenable by analysis techniques that assume the data to be normally distributed.  When setting
        mx.diff=FALSE , then Equation 4 is employed, calculating enrichment in an analogous way to classical
        GSEA which typically provides a bimodal distribution of GSVA enrichment scores for each gene.
    """

    if not single_sample and permutations:
        raise warn('permutations are not supported when not single_sample')

    key = (expression.hashable, method, gene_sets)

    if key in GSVA_CACHE:
        return GSVA_CACHE[key]

    if single_sample:
        assert isinstance(expression, Profile)
        joined = DataFrame(
            concat([expression.top.up, expression.top.down]),
        )
        joined.columns = ['condition']
        joined['control'] = 0
        joined.index = joined.index.astype(str)

        globalenv['expression'] = joined
        globalenv['expression_classes'] = Series(['case', 'control'])
    else:
        joined = expression.joined

        joined = DataFrame(joined)
        joined.index = joined.index.astype(str)

        nulls = joined.isnull().any(axis=0).reset_index(drop=True)
        if nulls.any():
            print(f'Following columns contain nulls and will be skipped: {list(joined.columns[nulls])}')
        joined = joined[joined.columns[~nulls]]

        globalenv['expression'] = joined
        globalenv['expression_classes'] = expression.classes.loc[~nulls.reset_index(drop=True)]

    mx_diff = 'T' if mx_diff else 'F'
    procedure = 'permutations' if single_sample else 'bayes'

    result = r(f"""
    result = gsva.with_probabilities(
        expression, expression_classes, {gene_sets}, '{procedure}',
        method = '{method}', mx.diff = {mx_diff}
        {', permutations = ' + str(permutations) if procedure == 'permutations' else ''}
    )
    rows = rownames(result)
    result
    """)
    result.index = r['rows']

    global process_specific_counter
    process_specific_counter += 1
    if process_specific_counter % 20 == 0:
        r('gc()')

    if _cache:
        GSVA_CACHE[key] = result
    return result


def create_gsva_scorer(
    gene_sets='c2.cp.kegg', id_type='entrez', grouping='by_substance',
    q_value_cutoff=0.1, na_action='fill_0', method='gsva', single_sample=False,
    permutations=None, mx_diff=True
):
    if single_sample and method == 'plage':
        warn('PLAGE is not suitable for single sample testing')

    importr('GSVA')
    importr('Biobase')
    gsea_base = importr('GSEABase')
    cwd = Path(__file__).parent
    r(f'source("{cwd}/gsva.R")')

    gmt_path = db.resolve(gene_sets, id_type)
    gene_sets_r = gsea_base.getGmt(gmt_path)

    input = Profile if single_sample else ExpressionWithControls

    def gsva_score(disease: input, compound: input):
        multiprocess_cache_manager.respawn_cache_if_needed()
        globalenv[gene_sets] = gene_sets_r

        disease_gene_sets = gsva(disease, gene_sets=gene_sets, method=method, single_sample=single_sample, permutations=permutations, mx_diff=mx_diff)

        disease_gene_sets.drop(disease_gene_sets[disease_gene_sets['fdr_q-val'] > q_value_cutoff].index, inplace=True)

        signature_gene_sets = gsva(compound, gene_sets=gene_sets, method=method, single_sample=single_sample, permutations=permutations, mx_diff=mx_diff, _cache=False)

        joined = combine_gsea_results(disease_gene_sets, signature_gene_sets, na_action)

        return joined.score.mean()

    gsva_score.__name__ = (
        method +
        f'_{permutations}' +
        f'_mx_diff:{mx_diff}' +
        ('_single_sample' if single_sample else '')
    )

    return scoring_function(gsva_score, input=input, grouping=grouping)
