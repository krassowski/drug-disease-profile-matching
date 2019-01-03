import sys
from abc import ABC, abstractmethod
from collections import UserDict
from typing import Tuple

from pandas import Series, concat, DataFrame, Index

from data_frames import AugmentedSeries, AugmentedDataFrame
from data_sources.drug_connectivity_map import SignaturesData, get_controls_for_signatures
from data_sources.tcga import ExpressionManager
from helpers import first


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


class ExpressionWithControls(ScoringInput):
    # TODO: make both Profile and ExpressionWithControls inherit from same base
    #  and have a limit() method (to trim the data to the most relevant genes only)

    cases: DataFrame
    controls: DataFrame
    joined: DataFrame
    classes: Series

    def __repr__(self):
        return f'<{self.__class__.__name__}: {len(self.cases.columns)} cases, {len(self.controls.columns)} controls>'

    @property
    def hashable(self):
        all_data = self.joined
        return tuple(all_data.index), tuple(all_data.columns), all_data.sum().sum()


class TCGAExpressionWithControls(ExpressionWithControls, ExpressionManager):

    @property
    def cases(self):
        data = self
        return data[data.columns[data.classes == 'tumor']]

    @property
    def controls(self):
        data = self
        return data[data.columns[data.classes == 'normal']]

    @property
    def joined(self):

        return self


class SignaturesWithControls(ExpressionWithControls, AugmentedDataFrame):

    @property
    def cases(self):
        return self

    @property
    def controls(self):
        controls = get_controls_for_signatures(tuple(self.columns), genes_to_keep=tuple(self.index))
        controls.columns = [c + '_control' for c in controls.columns]
        return controls

    @property
    def classes(self):
        return Series(['case'] * len(self.cases.columns) + ['normal'] * len(self.controls.columns))

    @property
    def joined(self):
        return concat([self.cases, self.controls], axis=1)


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


class SubstancesCollectionWithControls(SubstancesCollection, UserDict):
    """Signatures are grouped by the substance; each signature has it's set of controls"""

    def get_attr(self, key):
        values = self.values()
        f = getattr(first(values), key)
        for v in values:
            f2 = getattr(v, key)
            assert (f == f2).all()
        return f

    @property
    def genes(self):
        return self.get_attr('index')

    def groups_keys(self):
        return list(self.keys())

    @property
    def signature_ids(self):
        return {
            signature_id
            for signatures in self.values()
            for signature_id in signatures.columns
        }

    def drop_signatures(self, ids):
        return SubstancesCollectionWithControls({
            signature_ids: signatures.drop(columns=set(signatures.columns) & set(ids))
            for signature_ids, signatures in self.items()
        })

    @classmethod
    def from_signatures(cls, signatures: SignaturesData):
        data = {}
        for substance in set(signatures.classes()):
            signature_ids = signatures.members_of_class(substance)
            data[tuple(signature_ids)] = SignaturesWithControls(signatures[signature_ids])
        return cls(data)
