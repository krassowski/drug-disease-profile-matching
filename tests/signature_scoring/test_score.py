from pandas import Series
from pytest import approx

from signature_scoring.profile import Profile
from signature_scoring.profile import Signature
from signature_scoring.scoring_functions.generic_scorers import score_spearman


disease = Series({'BRCA1': 10, 'B': 1, 'T': -1, 'TP53': -10})
drug_1 = Series({'BRCA1': -10, 'T': -1, 'B': 1, 'TP53': +10})
drug_2 = Series({'BRCA1': -10, 'T': 1, 'B': -1, 'TP53': +10})
drug_3 = Series({'BRCA1': 2, 'TP53': -1, 'B': 1, 'T': -1})


def test_profile_ranks():

    disease_profile = Profile(disease)

    disease_ranks = {
        gene: i + 1
        for i, gene in enumerate(disease.sort_values().index)
    }
    #assert disease.rank().to_dict() == disease_ranks
    #assert disease_ranks == {'BRCA1': 4, 'B': 3, 'T': 2, 'TP53': 1}

    #assert disease_profile.top.ranks.to_dict() == disease_ranks
    #assert disease_profile.full.ranks.to_dict() == disease_ranks


def test_split():
    down, up = Signature(disease).split(10)
    assert list(down.index) == ['TP53', 'T']
    assert list(up.index) == ['BRCA1', 'B']


def test_profile_split():
    drug_profile = Profile(drug_1)

    drug_down, drug_up = Signature(drug_1).split(10)

    assert (drug_down == drug_profile.top.down).all()
    assert (drug_up == drug_profile.top.up).all()


def test_spearmanr():

    disease_profile = Profile(disease)
    drug_profile = Profile(drug_1)

    assert score_spearman(disease_profile, drug_profile) == approx(2)
    assert score_spearman(disease_profile, disease_profile) == approx(-2)
