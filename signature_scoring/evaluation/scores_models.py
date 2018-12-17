from typing import Dict

from dataclasses import dataclass
from pandas import Series, concat

from data_sources.drug_connectivity_map import AggregatedScores


@dataclass
class TopScores:
    all: Series
    indications: Series
    contraindications: Series
    controls: Series


@dataclass
class ScoresVector:
    """Collection of scores represented as vectors:
     - observed scores (scaled to [-1, +1]),
     - corresponding expected scores.
    """
    expected: list
    observed: Series

    def __init__(self, scores_map: Dict[int, AggregatedScores], limit_to=None):
        if limit_to:
            scores_map = {k: v for k, v in scores_map.items() if k in limit_to}

        perfect_result = []

        for weight, scores in scores_map.items():
            perfect_result += [weight] * len(scores)

        scores = [df.score for df in scores_map.values() if not df.empty]

        if not scores:
            print('Something went wrong, it might have been aggregation step')
            assert False

        all_results = concat(scores)
        result_scaled = (all_results - all_results.min()) / (all_results.max() - all_results.min())

        self.expected = perfect_result
        self.observed = -1 + 2 * result_scaled


@dataclass
class ProcessedScores:
    top: TopScores

    vector_overall: ScoresVector
    vector_controls: ScoresVector
    vector_contraindications: ScoresVector

    indications: Series
    contraindications: Series
    controls: Series

