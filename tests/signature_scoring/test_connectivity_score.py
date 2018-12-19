from pandas import Series, DataFrame

from signature_scoring import score_signatures
from signature_scoring.models import Profile
from signature_scoring.scoring_functions.connectivity_score import create_scorer


# query = dcm.from_perturbations(['vemurafenib'])
# query = query.head(n=25)
query = Series({
    b'5720': -0.16644155979156494,
    b'466': 0.46338778734207153,
    b'6009': 0.23561003804206848,
    b'2309': 0.14446471631526947,
    b'387': -0.5576758980751038,
    b'3553': -0.24347665905952454,
    b'427': 0.5373557806015015,
    b'5898': -0.5004006624221802,
    b'23365': -0.46615177392959595,
    b'6657': -0.22289252281188965,
    b'5054': 0.1859324425458908,
    b'3108': 1.1163060665130615,
    b'1950': -0.15650977194309235,
    b'351': 0.2328593134880066,
    b'4846': -0.3636510968208313,
    b'1452': 0.09088471531867981,
    b'4776': -0.2298683524131775,
    b'6908': 0.1434759646654129,
    b'672': 1.3136569261550903,
    b'5710': -0.08952534198760986,
    b'2115': -0.9770954847335815,
    b'7015': 0.0819830447435379,
    b'8726': 0.18349042534828186,
    b'2185': -0.2304331660270691,
    b'3315': 0.1746431589126587
})
query.name = 'CPC006_VCAP_24H:BRD-K56343971-001-02-3:10'


def test_connectivity_score():
    connectivity_score = create_scorer(negative=False)

    profile = Profile(query)

    score = connectivity_score(profile, profile)
    assert score > 0

    # negative=True: substances which might treat disease of given profile will get high scores
    reverse_connectivity_score = create_scorer(negative=True)

    score = reverse_connectivity_score(profile, profile)
    assert score < 0

    hypothetical_disease_profile = Profile(-query)
    ideal_drug_profile = Profile(query)

    score = reverse_connectivity_score(hypothetical_disease_profile, ideal_drug_profile)
    assert score > 0

    scores = score_signatures(connectivity_score, query, DataFrame(query), limit=None, processes=1)
    assert scores[query.name] > 0

    scores = score_signatures(reverse_connectivity_score, query, DataFrame(query), limit=None, processes=1)
    assert scores[query.name] < 0
