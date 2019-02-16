from numpy import sign, square
from pandas import Series

from helpers.inline import inline, compile_with_inline, inline_if_else

from ..models import Profile
from . import scoring_function


def kolmogorov_smirnov(ranks_disease, ranks_compound):

    a_list = []
    b_list = []

    n = len(ranks_disease)
    t = len(ranks_compound)

    for gene, j in ranks_compound.items():
        z = ranks_disease[gene] / n
        a_list.append((j + 1) / t - z)
        b_list.append(z - j / t)

    a = max(a_list)
    b = max(b_list)

    # assert (-1 <= a <= 1) and (-1 <= b <= 1)

    return a if a > b else -b


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


def kolmogorov_smirnov_plot(ranks_disease, ranks_compound):
    from seaborn import lineplot

    x_ecdf = []
    y_ecdf = []

    n = len(ranks_disease)
    t = len(ranks_compound)

    for gene, j in ranks_compound.to_dict().items():
        zi = ranks_disease[gene] / n
        x_ecdf.append(zi)
        y_ecdf.append(j / t)

    lineplot(range(t), x_ecdf)
    lineplot(range(t), y_ecdf)

    raise Exception


def difference(ks_up, ks_dn, factor):
    return (ks_up - ks_dn) * factor


def conditional_difference(ks_up, ks_dn, factor):
    if sign(ks_up) == sign(ks_dn):
        return 0
    else:
        return (ks_up - ks_dn) * factor


def max_up_or_down(ks_up, ks_dn, factor):
    return max(ks_up, ks_dn, key=abs) * factor


def create_scorer(negative, statistic=kolmogorov_smirnov, ranks_type='signature', compose_tags=conditional_difference):
    if negative:
        print('negative=True: substances which might treat disease of given profile will get high scores')
    else:
        print('negative=False: substances with similar effect on expression will get high scores')

    assert ranks_type in {'disease', 'signature'}

    """
    
    Cheng study:
    
    UpInDisease = a set of N up‐regulated features from disease genomic data
    DownInDisease = a set of N down‐regulated features from disease genomic data
    KSup = the KS score between UpInDisease and complete compound profile
    KSdown = the KS score between DownInDisease and complete compound profile
    IfKSup and KSdown have different signs then Connectivity score = KSup −  KSdown else Connectivity score = 0

    For this study, N is set to 500 in all metrics.
    """
    """
    As per CMap publication and more detailed:
    https://portals.broadinstitute.org/cmap/help_topics_linkified.jsp#how%20connectivity%20score%20is%20calculated
    """
    def connectivity_score(disease_profile: Profile, compound_profile: Profile):

        factor = inline_if_else(negative, -1, 1)
        if inline(ranks_type == 'disease'):
            ranks = disease_profile.full.ranks
            tags = {
                'up': compound_profile.top.up,
                'down': compound_profile.top.down
            }
        else:
            ranks = compound_profile.full.ranks
            tags = {
                'up': disease_profile.top.up,
                'down': disease_profile.top.down
            }

        compute_statistic = inline(statistic)

        ks = {}

        for tag_name, genes_in_tag in tags.items():
            ks[tag_name] = compute_statistic(ranks, genes_in_tag.rank(ascending=False))

        ks_dn, ks_up = ks['down'], ks['up']

        return inline(compose_tags(ks_up, ks_dn, factor))

    name = (
        'connectivity_score' +
        '_' + ranks_type +
        '_' + statistic.__name__ +
        '_' + compose_tags.__name__
    )
    return scoring_function(
        compile_with_inline(connectivity_score, name, locals(), globals())
    )
