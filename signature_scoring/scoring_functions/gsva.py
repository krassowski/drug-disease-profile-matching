from typing import Union

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


def gsva(expression: Union[ExpressionWithControls, Profile], gene_sets: str, method: str = 'gsva', single_sample=False):
    # TODO: p-value cutoff

    key = (expression.hashable, method, gene_sets)

    if key in GSVA_CACHE:
        return GSVA_CACHE[key]

    if isinstance(expression, Profile):
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

    # mx.diff
    """
    An important argument of the gsva() function is the flag mx.diff which is set to TRUE by default.
    
    Under this default setting, GSVA enrichment scores are calculated using Equation 5, and therefore, are
    more amenable by analysis techniques that assume the data to be normally distributed.  When setting
    mx.diff=FALSE , then Equation 4 is employed, calculating enrichment in an analogous way to classical
    GSEA which typically provides a bimodal distribution of GSVA enrichment scores for each gene.
    """

    result = r(f"""
    design = cbind(condition=expression_classes != 'normal')
    phenoData = AnnotatedDataFrame(data=as.data.frame(as.table(design)))
    row.names(phenoData) = colnames(expression)
    expression_set = ExpressionSet(assayData=data.matrix(expression), phenoData=phenoData)
    #  geneid=dimnames(expression)[[1]]
    result = gsva(exprs(expression_set), {gene_sets}, method='{method}', verbose=F, parallel.sz=1)

    gene_sets = rownames(result)
    samples = colnames(result)
    result
    """)
    rows = r['gene_sets']
    if single_sample:
        columns = r['samples']
        # result = DataFrame(data=result, index=rows, columns=joined.columns)
        result = DataFrame(data=result, index=rows, columns=columns)
        result['nes'] = result['condition']
        result['fdr_q-val'] = 0  # 1 - (result['condition'] - result['control']).abs()
    else:
        # result of the gsva is then used to create a table
        result = r("""
        library(limma)
        design = cbind(all=1, condition=expression_classes != 'normal')
        fit <- lmFit(result, design)
        fit <- eBayes(fit)
        allGeneSets <- topTable(fit, coef="condition", number=Inf, adjust="BH")
        # DeGeneSets <- topTable(fit, coef="condition", number=Inf, p.value=adjPvalueCutoff, adjust="BH")
        gene_sets = rownames(allGeneSets)
        allGeneSets
        """)
        rows = r['gene_sets']
        result = result.rename({'adj.P.Val': 'fdr_q-val', 'logFC': 'nes'}, axis=1)
        result.index = rows

    r('rm(expression_set, design, expression, result, gene_sets, samples); gc()')

    GSVA_CACHE[key] = result
    return result


def create_gsva_scorer(
    gene_sets='c2.cp.kegg', id_type='entrez', grouping='by_substance',
    q_value_cutoff=0.1, na_action='fill_0', method='gsva', single_sample=False
):

    importr('GSVA')
    importr('Biobase')
    gsea_base = importr('GSEABase')

    gmt_path = db.resolve(gene_sets, id_type)
    globalenv[gene_sets] = gsea_base.getGmt(gmt_path)

    input = Profile if single_sample else ExpressionWithControls

    def gsva_score(disease: input, compound: input):
        multiprocess_cache_manager.respawn_cache_if_needed()

        disease_gene_sets = gsva(disease, gene_sets=gene_sets, method=method, single_sample=single_sample)
        disease_gene_sets.drop(disease_gene_sets[disease_gene_sets['fdr_q-val'] > q_value_cutoff].index, inplace=True)

        signature_gene_sets = gsva(compound, gene_sets=gene_sets, method=method, single_sample=single_sample)

        joined = combine_gsea_results(disease_gene_sets, signature_gene_sets, na_action)
        return joined.score.mean()

    if single_sample:
        gsva_score.__name__ += '_single_sample'

    return scoring_function(gsva_score, input=input, grouping=grouping)
