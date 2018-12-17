import sys
from typing import Tuple

from pandas import Series, concat

from data_frames import AugmentedSeries


class Signature(AugmentedSeries):

    def split(self, limit: int, nlargest=Series.nlargest) -> Tuple[Series, Series]:
        up_regulated = self[self > 0]
        up_regulated = up_regulated[nlargest(up_regulated, limit).index].nlargest(limit)
        down_regulated = self[self < 0]
        down_regulated = down_regulated[nlargest(-down_regulated, limit).index].nsmallest(limit)
        return down_regulated, up_regulated


class ScoringInput:

    up: Series
    down: Series
    ranks: Series

    def __init__(self, signature: Signature, limit=None, nlargest=Series.nlargest):

        self.down, self.up = signature.split(limit or sys.maxsize, nlargest=nlargest)
        if limit:
            self.ranks = concat([self.up, self.down]).rank(ascending=False)
        else:
            self.ranks = signature.rank(ascending=False)

    @property
    def genes(self) -> set:
        return {*self.down.index, *self.up.index}


class Profile:

    top: ScoringInput
    full: ScoringInput

    def __init__(self, signature: Signature, limit=None, **kwargs):
        if not isinstance(signature, Signature):
            signature = Signature(signature)
        self.top = ScoringInput(signature, limit, **kwargs)
        self.full = ScoringInput(signature, **kwargs)
