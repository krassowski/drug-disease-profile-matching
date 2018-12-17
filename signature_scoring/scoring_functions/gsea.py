from functools import partial

from pandas import concat, Series

from multiprocess.cache_manager import multiprocess_cache_manager
from methods.gsea import GSEADesktop
from data_frames import AugmentedDataFrame
from models import ExpressionProfile

from ..profile import Profile

GSEA_CACHE = None
multiprocess_cache_manager.add_cache(globals(), 'GSEA_CACHE', 'dict')


class DummyExpressions(ExpressionProfile):

    @classmethod
    def from_differential(cls, differential, case_name='case', control_name='control'):
        """Signatures are already differential"""
        self = cls(differential)
        self.columns = [case_name]
        self['control'] = 0
        self.case_name = case_name
        self.control_name = control_name
        return self

    @property
    def classes(self):
        return Series([self.case_name, self.control_name])


def create_gsea_scorer(gsea_app=GSEADesktop, permutations=500, gene_sets='h.all', q_value_cutoff=0.1, na_action='fill_0', score='mean', normalization=True):
    """
    na_action: fill_0 or drop
    score: mean, max, sum
        # mean = mean improvement for the condition (balancing pros and cons)
        # max = best effect on a molecular pathway, though might lead to deterioration of certain pathways of the disease.
        # sum =
    """
    run_gsea_app = gsea_app().run

    def cached_gsea_run(gsea, gene_sets, q_value_cutoff, differential_profile, class_name, warn=False):

        profile_hash = hash(differential_profile)
        key = (gene_sets, q_value_cutoff, profile_hash)

        if key in GSEA_CACHE:
            results = GSEA_CACHE[key]
        else:
            if warn:
                print('Warning: not using cache (that\'s fine if it\'s the first run)')
            dummy_profile = DummyExpressions.from_differential(differential_profile, case_name=class_name)
            results = gsea(
                dummy_profile,
                outdir=f'.temp_gsea/{class_name}',
                name=f'{str(profile_hash).replace("-", "m")}'
            )
            GSEA_CACHE[key] = results

        return results[class_name], results['control']

    def gsea_score(disease_profile: Profile, compound_profile: Profile):
        # theoretically this could be replaced with the original data (e.g. tumour/normal for TCGA),
        # though this would be less consistent with how the signature profile is handled and require
        # re-calculation of the differential profile for each GSEA run.
        multiprocess_cache_manager.respawn_cache_if_needed()

        disease_differential = AugmentedDataFrame(
            concat([disease_profile.top.up, disease_profile.top.down])
        )

        gsea = partial(
            run_gsea_app,
            gene_sets=gene_sets, id_type='entrez',
            permutations=permutations, permutation_type='Gene_set',
            metric='Diff_of_Classes',
            normalization='meandiv' if normalization else None
        )

        disease_gene_sets_up, disease_gene_sets_dn = cached_gsea_run(
            gsea, gene_sets, q_value_cutoff, disease_differential, class_name='disease', warn=True
        )

        signature = AugmentedDataFrame(
            concat([compound_profile.top.up, compound_profile.top.down])
        )

        signature_gene_sets_up, signature_gene_sets_dn = cached_gsea_run(
            gsea, gene_sets, q_value_cutoff, signature, class_name='signature'
        )

        results = [disease_gene_sets_up, disease_gene_sets_dn]

        if q_value_cutoff:
            for result in results:
                result.drop(result[result['fdr_q-val'] > q_value_cutoff].index, inplace=True)

        disease_gene_sets = concat([disease_gene_sets_up, disease_gene_sets_dn])
        signature_gene_sets = concat([signature_gene_sets_up, signature_gene_sets_dn])

        joined = disease_gene_sets.merge(
            signature_gene_sets, suffixes=['_disease', '_signature'],
            left_on=disease_gene_sets.index, right_on=signature_gene_sets.index,
            how=('left' if na_action == 'fill_0' else 'inner')
        )
        joined = joined.set_index('key_0')
        #from utilities_namespace import show_table
        #show_table(joined)

        if na_action == 'fill_0':
            joined = joined.fillna(0)

        if normalization:
            joined['score'] = (
                1
                -
                (joined.nes_disease + joined.nes_signature) / joined.nes_disease
                *
                (1 - joined['fdr_q-val_disease']) * (1 - joined['fdr_q-val_disease'])
            )
        else:
            joined['score'] = (
                1
                -
                (joined.es_disease + joined.es_signature) / joined.es_disease
                *
                (1 - joined['fdr_q-val_disease']) * (1 - joined['fdr_q-val_disease'])
            )

        return getattr(joined.score, score)()

    gsea_score.__name__ = f'gsea_{permutations}_{gene_sets}_{score}_{q_value_cutoff}_{na_action}_{normalization}'

    return gsea_score

