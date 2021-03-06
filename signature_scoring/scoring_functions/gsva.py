from pathlib import Path
from typing import Union
from warnings import warn
from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE

from numpy import nan
from pandas import read_csv
from pandas import DataFrame, concat, Series
from pandas.errors import EmptyDataError
from rpy2.robjects import r
from rpy2.robjects.packages import importr

from helpers.temp import create_tmp_dir
from data_sources.molecular_signatures_db import MolecularSignaturesDatabase
from enhanced_multiprocessing.cache_manager import multiprocess_cache_manager

from ..models import Profile
from ..models.with_controls import ExpressionWithControls
from . import scoring_function
from .gsea import combine_gsea_results


db = MolecularSignaturesDatabase()
GSVA_CACHE = None
multiprocess_cache_manager.add_cache(globals(), 'GSVA_CACHE', 'dict')


gsva_tmp_dir = create_tmp_dir('gsva')
vanilla_R = ['R', '--vanilla', '--quiet']


def gsva(
    expression: Union[ExpressionWithControls, Profile], gene_sets_path: str, method: str = 'gsva',
    single_sample=False, permutations=1000, mx_diff=True, cores=1, _cache=True, limit_to_gene_sets=False,
    verbose=False
):
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

    key = (expression.hashable, method, gene_sets_path)

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

        expression_classes = Series(['case', 'control'])
        expression = joined

    else:
        joined = expression.joined

        joined = DataFrame(joined)
        joined.index = joined.index.astype(str)

        nulls = joined.isnull().any(axis=0).reset_index(drop=True)
        if nulls.any():
            print(f'Following columns contain nulls and will be skipped: {list(joined.columns[nulls])}')
        joined = joined[joined.columns[~nulls]]

        expression_classes = expression.classes.loc[~nulls.reset_index(drop=True)]
        expression = joined

    mx_diff = 'T' if mx_diff else 'F'
    procedure = 'gene_permutation' if single_sample else 'bayes'
    cwd = Path(__file__).parent

    with NamedTemporaryFile(mode='w', prefix=gsva_tmp_dir) as f_expression, NamedTemporaryFile(prefix=gsva_tmp_dir) as f_result:
        expression.to_csv(f_expression)
        script = f"""
        source("{cwd}/gsva.R")
        expression = read.csv('{f_expression.name}', row.names=1)
        expression_classes = c{tuple(expression_classes)}
        gene_sets = readRDS('{gene_sets_path}')

        result = gsva.with_probabilities(
            expression, expression_classes, gene_sets, '{procedure}',
            method = '{method}', mx.diff={mx_diff}, include_control=F, cores={cores},
            limit_to_gene_sets={'c' + str(tuple(limit_to_gene_sets)) if limit_to_gene_sets is not False else 'F'}, progress=F
            {', permutations = ' + str(permutations) if procedure == 'permutations' else ''}
        )
        write.csv(result, '{f_result.name}')
        """

        process = Popen(vanilla_R, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        r = process._stdin_write(script.encode())
        if verbose:
            from helpers.streams import handle_streams
            from signature_scoring.evaluation import display

            handlers = {'out': display, 'err': warn}
            handle_streams(process, handlers)
        else:
            process.wait()
        try:
            result = read_csv(f_result.name, index_col=0)
        except EmptyDataError:
            result = DataFrame()

    if _cache:
        GSVA_CACHE[key] = result
    return result


def create_gsva_scorer(
    gene_sets='c2.cp.kegg', id_type='entrez', grouping='by_substance',
    q_value_cutoff=0.1, na_action='fill_0', method='gsva', single_sample=False,
    permutations=None, mx_diff=True, custom_multiprocessing=False
):

    # as long as dummy controls are not included (include_control=F) there is no problem, otherwise:
    # if single_sample and method == 'plage':
    #    warn('PLAGE is not suitable for single sample testing')

    gmt_path = db.resolve(gene_sets, id_type)
    gsea_base = importr('GSEABase')
    gene_sets_r = gsea_base.getGmt(gmt_path)

    # transform to named list from GeneSetCollection class object
    gene_sets_r = gsea_base.geneIds(gene_sets_r)

    gene_sets_file = NamedTemporaryFile(delete=False, prefix=gsva_tmp_dir)
    r['saveRDS'](gene_sets_r, file=gene_sets_file.name)

    input = Profile if single_sample else ExpressionWithControls

    def gsva_score(disease: input, compound: input, cores=1):
        if not custom_multiprocessing:
            multiprocess_cache_manager.respawn_cache_if_needed()

        disease_gene_sets = gsva(
            disease, gene_sets_path=gene_sets_file.name, method=method, single_sample=single_sample,
            permutations=permutations, mx_diff=mx_diff, cores=cores
        )

        disease_gene_sets.drop(disease_gene_sets[disease_gene_sets['fdr_q-val'] > q_value_cutoff].index, inplace=True)

        assert len(disease_gene_sets.index)
        signature_gene_sets = gsva(
            compound, gene_sets_path=gene_sets_file.name, method=method, single_sample=single_sample,
            permutations=permutations, mx_diff=mx_diff, _cache=False, cores=cores,
            limit_to_gene_sets=list(disease_gene_sets.index)
        )
        if signature_gene_sets.empty:
            return nan

        joined = combine_gsea_results(disease_gene_sets, signature_gene_sets, na_action)
        assert not (set(signature_gene_sets.index) - set(disease_gene_sets.index))
        return joined.score.mean()

    gsva_score.__name__ = (
        method +
        f'_{permutations}' +
        f'_mx_diff:{mx_diff}' +
        ('_single_sample' if single_sample else '')
    )

    return scoring_function(gsva_score, input=input, grouping=grouping, custom_multiprocessing=custom_multiprocessing)
