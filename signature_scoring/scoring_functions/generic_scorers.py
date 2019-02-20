from pandas import concat
from scipy.stats import spearmanr
from scipy import spatial

from ..models import Profile, ScoringData
from . import scoring_function


@scoring_function
def score_spearman(disease_profile: Profile, compound_profile: Profile):

    down_ranks = compound_profile.top.down.index
    up_ranks = compound_profile.top.up.index

    corresponding_ranks = disease_profile.top.ranks[up_ranks]
    signature_ranks = compound_profile.top.up.rank()
    s_up = spearmanr(signature_ranks, corresponding_ranks).correlation

    corresponding_ranks = disease_profile.top.ranks[down_ranks]
    signature_ranks = compound_profile.top.down.rank()
    s_dn = spearmanr(signature_ranks, corresponding_ranks).correlation

    return s_dn + s_up


@scoring_function
def score_spearman_max(disease_profile: Profile, compound_profile: Profile):

    down_ranks = compound_profile.top.down.index
    up_ranks = compound_profile.top.up.index

    corresponding_ranks = disease_profile.top.ranks[up_ranks]
    signature_ranks = compound_profile.top.up.rank()
    s_up = spearmanr(signature_ranks, corresponding_ranks).correlation

    corresponding_ranks = disease_profile.top.ranks[down_ranks]
    signature_ranks = compound_profile.top.down.rank()
    s_dn = spearmanr(signature_ranks, corresponding_ranks).correlation

    return max([s_dn, s_up])


def changed_subsets(disease: ScoringData, compound: ScoringData):
    changed_by_compound_series = concat([compound.down, compound.up])
    changed_by_compound = set(changed_by_compound_series.index)

    x_up_in_disease = changed_by_compound.intersection(disease.up.index)
    x_down_in_disease = changed_by_compound.intersection(disease.down.index)

    return changed_by_compound_series, x_down_in_disease, x_up_in_disease


@scoring_function
def x_sum(disease_profile: Profile, compound_profile: Profile):
    """Algorithm:

    XSum: The eXtreme Sum score is calculated as follows.

    UpInDisease and DownInDisease are defined as above

    ChangedByCompound = top N up‐regulated and N down‐regulated features by fold change values between compound treated samples and control samples

    XUpInDisease = UpInDisease ∩ ChangedByCompound
    XDownInDisease = DownInDisease ∩ ChangedByCompound

    sum(XUpInDisease) = sum of compound gene expression fold change values in the set XUpInDisease
    sum(XDownInDisease) = sum of compound gene expression fold change values in the set XDownInDisease
    XSum = sum(XUpInDisease) − sum(XDownInDisease)

    Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4278345/
    """

    changed_by_compound, x_down_in_disease, x_up_in_disease = changed_subsets(
        disease_profile.top, compound_profile.top
    )

    return changed_by_compound[x_down_in_disease].sum() - changed_by_compound[x_up_in_disease].sum()


@scoring_function
def x_sum_max(disease_profile: Profile, compound_profile: Profile):

    changed_by_compound, x_down_in_disease, x_up_in_disease = changed_subsets(
        disease_profile.top, compound_profile.top
    )

    return max([
        changed_by_compound[x_down_in_disease].sum(),
        - changed_by_compound[x_up_in_disease].sum()
    ])


@scoring_function
def x_product(disease_profile: Profile, compound_profile: Profile):

    changed_by_compound, x_down_in_disease, x_up_in_disease = changed_subsets(
        disease_profile.top, compound_profile.top
    )

    disease = disease_profile.top

    return (
        - (changed_by_compound[x_down_in_disease] * disease.down[x_down_in_disease]).sum()
        - (changed_by_compound[x_up_in_disease] * disease.up[x_up_in_disease]).sum()
    )


@scoring_function
def x_cos(disease_profile: Profile, compound_profile: Profile):
    changed_by_compound, x_down_in_disease, x_up_in_disease = changed_subsets(
        disease_profile.top, compound_profile.top
    )

    full_disease = concat([disease_profile.top.up, disease_profile.top.down])
    full_compound = concat([compound_profile.top.up, compound_profile.top.down])

    changed_in_disease = set(full_disease.index)

    x_set = changed_in_disease.intersection(changed_by_compound.index)

    disease = full_disease[full_disease.index[full_disease.index.isin(x_set)]]
    drug = full_compound.loc[disease.index]

    return spatial.distance.cosine(drug.values, disease.values)


@scoring_function
def x_product_max(disease_profile: Profile, compound_profile: Profile):

    changed_by_compound, x_down_in_disease, x_up_in_disease = changed_subsets(
        disease_profile.top, compound_profile.top
    )

    disease = disease_profile.top

    return max([
        - (changed_by_compound[x_down_in_disease] * disease.down[x_down_in_disease]).sum(),
        - (changed_by_compound[x_up_in_disease] * disease.up[x_up_in_disease]).sum()
    ])
