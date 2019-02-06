import sys
from abc import ABC, abstractmethod
from typing import Tuple

from pandas import Series, concat, Index

from data_frames import AugmentedSeries, AugmentedDataFrame


class Signature(AugmentedSeries):

    def split(self, limit: int, nlargest=Series.nlargest) -> Tuple[Series, Series]:
        up_regulated = self[self > 0]
        up_regulated = up_regulated[nlargest(up_regulated, limit).index].nlargest(limit)
        down_regulated = self[self < 0]
        down_regulated = down_regulated[nlargest(-down_regulated, limit).index].nsmallest(limit)
        return down_regulated, up_regulated


class ScoringData:

    up: Series
    down: Series
    ranks: Series

    def __init__(self, signature: Signature, limit=None, nlargest=Series.nlargest):

        self.down, self.up = signature.split(limit or sys.maxsize, nlargest=nlargest)
        if limit:
            self.ranks = concat([self.up, self.down]).rank(ascending=False)
        else:
            self.ranks = signature.rank(ascending=False)
        self._hashable = tuple(tuple(component.items()) for component in [self.up, self.down, self.ranks])

    @property
    def genes(self) -> set:
        return {*self.down.index, *self.up.index}

    @property
    def hashable(self):
        return self._hashable


class ScoringInput(ABC):

    @property
    @abstractmethod
    def hashable(self):
        pass


class Profile(ScoringInput):

    top: ScoringData
    full: ScoringData

    def __init__(self, signature: Signature, limit=None, **kwargs):
        if not isinstance(signature, Signature):
            signature = Signature(signature)
        self.top = ScoringData(signature, limit, **kwargs)
        self.full = ScoringData(signature, **kwargs)

    @property
    def hashable(self):
        return self.top.hashable, self.full.hashable


class SignaturesGrouping(ABC):

    @property
    @abstractmethod
    def signature_ids(self) -> set:
        """Return all signatures from all the groups"""

    @abstractmethod
    def groups_keys(self) -> list:
        """Return names of all groups"""

    @property
    @abstractmethod
    def genes(self) -> Index:
        pass

    @abstractmethod
    def drop_signatures(self, ids):
        pass


class SignaturesCollection(SignaturesGrouping, AugmentedDataFrame):
    """Each signature has it's own group"""

    def groups_keys(self):
        return list(self.columns)

    @property
    def signature_ids(self):
        return set(self.columns)

    @property
    def genes(self):
        return self.index

    def drop_signatures(self, ids):
        return self.drop(columns=ids)


class SubstancesCollection(SignaturesGrouping):
    """Signatures are grouped by the substance"""
