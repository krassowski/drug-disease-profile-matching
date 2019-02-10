from collections import UserDict
from copy import copy
from typing import TextIO

from pandas import DataFrame, Series, concat, np

from data_frames import AugmentedDataFrame
from data_sources.drug_connectivity_map import get_controls_for_signatures, SignaturesData
from data_sources.tcga import ExpressionManager
from helpers import first
from models import ExpressionProfile

from . import ScoringInput, SubstancesCollection


class ExpressionWithControls(ScoringInput, ExpressionProfile):
    # TODO: have a limit() method (to trim the data to the most relevant genes only)

    cases: DataFrame
    controls: DataFrame
    joined: DataFrame
    classes: Series
    case_name: str
    control_name: str

    def __repr__(self):
        return f'<{self.__class__.__name__}: {len(self.cases.columns)} cases, {len(self.controls.columns)} controls>'

    @property
    def genes(self) -> Series:
        return self.index

    @property
    def hashable(self):
        all_data = self.joined
        return tuple(all_data.index), tuple(all_data.columns), all_data.sum().sum()

    def to_cls(self, f):
        classes = self.classes
        f.write(f'{len(classes)} {len(set(classes))} 1\n')
        classes_set = []
        for class_ in classes:
            if class_ not in classes_set:
                classes_set.append(class_)
        f.write(f'# {" ".join(classes_set)}\n')
        f.write(' '.join(classes))

    def to_gct(self, f: TextIO, tabular_writer='to_txt'):
        f.write('#1.2\n')
        expression_data = self.joined
        assert expression_data.notnull().all().all()
        f.write(f'{len(expression_data)}\t{len(expression_data.columns)}\n')
        getattr(self, tabular_writer)(f, expression_data)

    def to_txt(self, f: TextIO, expression_data=None):
        if expression_data is None:
            expression_data = self.joined
        columns = expression_data.columns
        expression_data['Description'] = 'na'
        expression_data = expression_data[['Description', *columns]]
        if type(expression_data.index[0]) is bytes:
            expression_data.index = [b.decode('utf-8') for b in expression_data.index]
        expression_data.index = expression_data.index.astype(int)
        expression_data.index.name = 'gene'
        header = '\t'.join([expression_data.index.name, *expression_data.columns]) + '\n'
        f.write(header)

        # 13.1 s
        np.savetxt(
            f,
            expression_data.reset_index().values,
            delimiter='\t',
            # entrez_id (integer), 'na' (not a description, str), *data (floats)
            fmt='%d\t%s' + ('\t%f' * (len(expression_data.columns) - 1))
        )

        # 19.6 s
        # expression_data = expression_data.astype(str)
        # for l in expression_data.itertuples():
        #    f.write('\t'.join(l) + '\n')

        # 19.3 s
        # lines = []
        # for l in expression_data.itertuples():
        #     lines.append('\t'.join(l))
        # f.write('\n'.join(lines))

        # 21.6 s
        # for l in expression_data.index.str.cat(expression_data, sep='\t'):
        #     f.write(l + '\n')

        # 21.1 s
        #  for l in zip(*[expression_data.index, *[expression_data[column] for column in expression_data.columns]]):
        #    f.write('\t'.join(l) + '\n')

        # 36.7 s
        # expression_data.to_csv(f.name, sep='\t', mode='a')

    def to_txt_naive(self, f: TextIO, expression_data=None):
        if expression_data is None:
            expression_data = self.joined

        columns = expression_data.columns
        expression_data['Description'] = 'na'
        expression_data = expression_data[['Description', *columns]]
        if type(expression_data.index[0]) is bytes:
            expression_data.index = [b.decode('utf-8') for b in expression_data.index]
        expression_data.index.name = 'gene'
        expression_data.to_csv(f, sep='\t')


class DummyExpressionsWithControls(ExpressionWithControls):

    @classmethod
    def from_differential(cls, differential, case_name='case', control_name='control'):
        """Signatures are already differential"""
        self = cls(differential)
        self.columns = [case_name]
        self['control'] = 0
        self.case_name = case_name
        self.control_name = control_name
        return self

    @property
    def classes(self):
        return Series([self.case_name, self.control_name])

    @property
    def joined(self):
        # the below casting might not be needed
        return DummyExpressionsWithControls(copy(self))


class TCGAExpressionWithControls(ExpressionWithControls, ExpressionManager):

    case_name = 'tumor'
    control_name = 'normal'

    @property
    def cases(self):
        data = self
        return data[data.columns[data.classes == self.case_name]]

    @property
    def controls(self):
        data = self
        return data[data.columns[data.classes == self.control_name]]

    @property
    def joined(self):
        return TCGAExpressionWithControls(copy(self))


class SignaturesWithControls(ExpressionWithControls, AugmentedDataFrame):

    case_name = 'case'
    control_name = 'normal'
    _metadata = ['_controls_cache']

    @property
    def cases(self):
        return self

    @property
    def controls(self):
        if hasattr(self, '_controls_cache') and self._controls_cache is not None:
            return self._controls_cache
        controls = get_controls_for_signatures(tuple(self.columns), genes_to_keep=tuple(self.index))
        controls.columns = [c + '_control' for c in controls.columns]
        self._controls_cache = controls
        return controls

    @property
    def classes(self):
        return Series(
            [self.case_name] * len(self.cases.columns)
            +
            [self.control_name] * len(self.controls.columns)
        )

    @property
    def joined(self):
        assert self.controls.notnull().all().all()
        return concat([self.cases, self.controls], axis=1)

    def __repr__(self):
        return DataFrame.__repr__(self)


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
