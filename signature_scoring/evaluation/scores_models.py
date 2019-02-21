from typing import List
from dataclasses import dataclass

from pandas import Series, DataFrame


Group = str  # in indications, controls, contraindications


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

    def __init__(self, scores: DataFrame, limit_to: List[Group] = None, rescale=True):
        if limit_to:
            scores = scores[scores.group.isin(limit_to)]

        if scores.empty:
            print('Something went wrong, it might have been aggregation step')
            assert False

        perfect_result = scores.expected_score
        observed = scores.score

        if rescale:
            observed_scaled = (observed - observed.min()) / (observed.max() - observed.min())
            observed = -1 + 2 * observed_scaled

        self.expected = perfect_result
        self.observed = observed


@dataclass
class ProcessedScores:
    top: TopScores

    vector_overall: ScoresVector
    vector_controls: ScoresVector
    vector_contraindications: ScoresVector

    vector_indications_over_non_indications: ScoresVector

    vector_overall_binary: ScoresVector
    vector_controls_binary: ScoresVector
    vector_contraindications_binary: ScoresVector

    indications: Series
    contraindications: Series
    controls: Series
    unassigned: Series
