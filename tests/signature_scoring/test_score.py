from pandas import Series
from pytest import approx

from signature_scoring.models import Profile
from helpers.mathtools import split_to_pos_and_neg
from helpers.cache import hash_series
from signature_scoring.models import Signature
from signature_scoring.scoring_functions.generic_scorers import score_spearman
import signature_scoring.scoring_functions.connectivity_score as connectivity


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


def test_zero_split():

    pos, neg = split_to_pos_and_neg(disease.values)
    assert list(pos) == [0, 1]
    assert list(neg) == [2, 3]

    pos, neg = split_to_pos_and_neg(Series({'A': 1, 'B': 0, 'C': 1}).values)
    assert list(pos) == [0, 2]
    assert list(neg) == []


def test_split():
    down, up = Signature(disease).split(10)
    assert list(down.index) == ['TP53', 'T']
    assert list(up.index) == ['BRCA1', 'B']

    down, up = Signature(Series(
        dict(zip(
            'ABCDEFGHI',
            range(1, 10)
        ))
    )).split(5)

    assert list(down.index) == []
    assert list(up.index) == ['I', 'H', 'G', 'F', 'E']


def test_hash():
    # identical objects have same hash
    assert hash_series(drug_1) == hash_series(drug_1)
    # different values make different hashes
    assert hash_series(drug_1) != hash_series(drug_2)
    # different indexes make different hashes
    drug_1_reindexed = drug_1.reindex(['B', 'BRCA1', 'TP53', 'T'])
    assert hash_series(drug_1) != hash_series(drug_1_reindexed)


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


def test_cramer():
    assert connectivity.cramér_von_mises(disease, drug_1) == 50
    assert connectivity.cramér_von_mises(disease, disease) == 0
