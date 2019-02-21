from copy import copy
from functools import partial
from typing import Dict

from numpy import sign, square
from pandas import Series, concat
from rpy2.robjects import r

from helpers.inline import inline, compile_with_inline, inline_if_else

from ..models import Profile
from . import scoring_function


def generalized_kolmogorov_smirnov(instance: Series, tag_list: Series, p: float = 1):
    """Random walk / running sums statistic of GSEA, generalised KS,

    with p=0 kind-of equivalent to normal KS (see proof in GSEA 2005 paper supplement)
    """
    ranked_list = instance.rank(ascending=False) + 1

    n = len(ranked_list)
    t = len(tag_list)

    decrement = 1 / (n - t)

    hit_denominator = ranked_list[tag_list.index].pow(p).sum()

    running_sum_statistic_hits = Series(data=0, index=instance.index)
    running_sum_statistic_hits[tag_list.index] = ranked_list[tag_list.index].pow(p) / hit_denominator
    running_sum_statistic_hits = running_sum_statistic_hits.cumsum()

    running_sum_statistic_misses = Series(data=decrement, index=instance.index)
    running_sum_statistic_misses[tag_list.index] = 0
    running_sum_statistic_misses = running_sum_statistic_misses.cumsum()

    return - max(running_sum_statistic_hits - running_sum_statistic_misses, key=abs)


def create_generalized_kolmogorov_smirnov(p: float = 1):
    compute_ks = partial(generalized_kolmogorov_smirnov, p=p)
    compute_ks.__name__ = f'generalized_kolmogorov_smirnov_{p}'
    return compute_ks


def kolmogorov_smirnov(
    instance: Series, tag_list: Series,
    zero_based_j=True, visualize=False, check_order=True, divide_by_n=True
):
    """See:
    - https://portals.broadinstitute.org/cmap/help_topics_linkified.jsp#how%20connectivity%20score%20is%20calculated
    - https://github.com/Bin-Chen-Lab/RGES/blob/af7fa59649dbc5dfa21b59f5446ea743ba236438/RGES_IC50_CMap_data.R

    for reference. Other loosely related resources:
    - https://github.com/cmap/psp/blob/master/broadinstitute_psp/sip/sip.py
    - https://github.com/cmap/merino/blob/4d5df231464528dc6719f74ab6fe86cc5a1f52a3/compute_wtcs.py
    - https://github.com/ewouddt/CSFA/blob/master/vignettes/CSFA.pdf

    """

    if check_order:
        assert (instance.values == instance.sort_values(ascending=False).values).all()

    n = len(instance)
    t = len(tag_list)

    misses_denominator = inline_if_else(divide_by_n, n, (n - t))

    tag_index_sorted_by_instance = instance.index[instance.index.isin(tag_list.index)]

    assert len(tag_list.index)

    if len(tag_index_sorted_by_instance) < 1:
        return 0

    # 1-based j over t
    running_sum_statistic_hits = Series(data=1 / t, index=tag_index_sorted_by_instance)
    running_sum_statistic_hits_cum = running_sum_statistic_hits.cumsum()

    # according to the paper, it should be just divided by n,
    # not by (n - t), but this give worse correlation with CMap 2.0 scores!
    # (and n - t is the GSEA way, and they say that their approach is based on GSEA KS, so this sounds ok)
    running_sum_statistic_misses_cum = (
        instance.rank(ascending=False) + 1
    )[tag_index_sorted_by_instance] / misses_denominator

    if visualize:
        kolmogorov_smirnov_plot(running_sum_statistic_misses_cum, running_sum_statistic_hits_cum)

    if zero_based_j:
        a = (running_sum_statistic_hits_cum - running_sum_statistic_misses_cum).max()

        # 0-based j over t
        running_sum_statistic_hits.iloc[0] = 0
        running_sum_statistic_hits_cum = running_sum_statistic_hits.cumsum()

        b = (running_sum_statistic_misses_cum - running_sum_statistic_hits_cum).max()
        return -(a if a > b else -b)
    else:
        return - max(running_sum_statistic_hits_cum - running_sum_statistic_misses_cum, key=abs)


def create_kolmogorov_smirnov(proper_ks=True, visualize=False, check_order=True, denominator='n'):

    assert denominator in ['n-t', 'n']

    compute_ks = partial(
        kolmogorov_smirnov,
        zero_based_j=proper_ks,
        visualize=visualize,
        check_order=check_order,
        divide_by_n=denominator == 'n'
    )

    compute_ks.__name__ = (
        f'kolmogorov_smirnov'
        f'_{"safe" if check_order else "fast"}'
        f'_{"properKS" if proper_ks else "KSbased"}'
        f'_{denominator.replace("-", "")}'
    )
    return compute_ks


def kolmogorov_smirnov_plot(running_sum_statistic_misses, running_sum_statistic_hits):
    from seaborn import lineplot
    from matplotlib import pyplot as plt
    plt.figure()
    lineplot(range(t), running_sum_statistic_misses, drawstyle='steps-pre', label='V(j)/n')
    lineplot(range(t), running_sum_statistic_hits, drawstyle='steps-pre', label='j/t')


def summation(ranks_disease, ranks_compound):
    running_sum = 0

    best = 0
    best_abs = 0

    n = len(ranks_disease)
    t = len(ranks_compound)

    for gene, j in ranks_compound.items():
        running_sum += j / t - ranks_disease[gene] / n
        if abs(running_sum) > best_abs:
            best = running_sum
            best_abs = abs(running_sum)

    return best


def cramér_von_mises(ranks_disease, ranks_compound: Series):

    n = len(ranks_disease)
    t = len(ranks_compound)

    # benchmarked on real data

    # naive approach: 229 ms
    # return sum(
    #     (j / t - ranks_disease[gene] / n) ** 2
    #     for gene, j in ranks_compound.items()
    # )

    # vectorized: 2.5 ms
    # return ((ranks_compound / t - ranks_disease[ranks_compound.index] / n) ** 2).sum()

    # vectorized + numpy: 2.26 ms
    return square(ranks_compound / t - ranks_disease[ranks_compound.index] / n).sum()


def difference(ks_up, ks_dn, factor):
    return (ks_up - ks_dn) * factor


def conditional_difference(ks_up, ks_dn, factor):
    if sign(ks_up) == sign(ks_dn):
        return 0
    else:
        return (ks_up - ks_dn) * factor


def max_up_or_down(ks_up, ks_dn, factor):
    return max(ks_up, ks_dn, key=abs) * factor


def create_scorer(
    negative,
    statistic=kolmogorov_smirnov,
    ranks_type='signature',
    compose_tags=conditional_difference,
    force_custom_tags: Dict[str, Series] = None
):
    if negative:
        print('negative=True: substances which might treat disease of given profile will get high scores')
    else:
        print('negative=False: substances with similar effect on expression will get high scores')

    assert ranks_type in {'disease', 'signature'}

    """
    As defined in Cheng, 2014 study and based on CMap publication, using the detailed description from:
    https://portals.broadinstitute.org/cmap/help_topics_linkified.jsp#how%20connectivity%20score%20is%20calculated
    """
    compute_statistic = statistic

    def connectivity_score(disease_profile: Profile, compound_profile: Profile):

        factor = inline_if_else(negative, -1, 1)
        if inline(ranks_type == 'disease'):
            ranks = disease_profile.full.ranks
            tags = {
                'up': compound_profile.full.up,
                'down': compound_profile.full.down
            }
        else:
            ranks = compound_profile.full.ranks
            tags = {
                'up': disease_profile.top.up,
                'down': disease_profile.top.down
            }

        if inline(force_custom_tags):
            # ignore tags inferred from the disease, use provided ones
            ranks = compound_profile.full.up
            tags = force_custom_tags
        else:
            pass

        # compute_statistic = inline(statistic)
        ranks = ranks.sort_values(ascending=False)

        ks = {}

        for tag_name, genes_in_tag in tags.items():
            ks[tag_name] = compute_statistic(
                ranks,
                genes_in_tag
            )

        ks_dn, ks_up = ks['down'], ks['up']

        return inline(compose_tags(ks_up, ks_dn, factor))

    name = (
        'connectivity_score' +
        '_' + ranks_type +
        '_' + statistic.__name__ +
        '_' + compose_tags.__name__
    )
    return scoring_function(
        compile_with_inline(connectivity_score, name, copy(locals()), {'force_custom_tags': force_custom_tags, **globals(), **locals()})
    )


r("""
connectivity_from_pharmaco_gx <- function(disease, drug, ...) {
    # full drug (compound) profile against top/down of the disease expression
    drug_data <- cbind(drug[, 2])
    rownames(drug_data) <- drug[, 1]
    
    rownames(disease) <- disease[, 1]
    disease_data <- disease[-1]

    PharmacoGx::connectivityScore(drug_data, disease_data, method='gsea', nperm=100)
}
""")


@scoring_function
def pharmaco_gx_connectivity_score(disease, drug):
    """PharmacoGx is distributed under the terms of the Artistic License 2.0"""
    drug_full = concat([drug.full.up, drug.full.down])
    disease_up_and_down = concat([disease.top.up, disease.top.down])

    for profile in [drug_full, disease_up_and_down]:
        profile.index = profile.index.astype(str)

    drug_full = drug_full.reset_index()
    disease_up_and_down = disease_up_and_down.reset_index()

    result = r['connectivity_from_pharmaco_gx'](drug_full, disease_up_and_down)
    result_data = dict(result.items())

    if 'es' in result_data:
        return result_data['es'] * -1
    else:
        return result_data['score'] * -1
