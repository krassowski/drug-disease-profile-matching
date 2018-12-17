from pronto import Ontology
from pandas import read_json, read_table, Series, DataFrame
from data_frames import MyDataFrame
from collections import Counter
from helpers.cache import cached_property
from gui_helpers import HorizontalNamespace
from helpers.mathtools import identity
from itertools import permutations
from functools import lru_cache
from warnings import warn
from typing import Dict
from copy import copy
from numpy import nan


class SRAAccessionsMap:
    
    def __init__(self, path='data/sra/SRA_Accessions.tab.xz'):
        self.path = path
        self.head = read_table(path, nrows=10)#, skiprows=25000)
        self.columns = self.head.columns
        
        # create lazy mapping options
        for a, b in permutations(self.columns, 2):
            lazy_map = self.lazy_map(a, b)
            name = (a + '_to_' + b).lower()
            setattr(self, name, lazy_map)

    def lazy_map(self, key, value):
        def closure(coerce=False):
            return self.map(key, value, coerce)
        return closure
    
    @lru_cache()
    def map(self, key, value, coerce=False) -> Dict:
        minimal_df = read_table(
            self.path, usecols=[key, value],
            na_values=['-'],
        ).drop_duplicates()
        
        skipped = minimal_df[minimal_df.isna().any(axis='columns')]
        if skipped.any().any():
            warn('Skipping NA values:')
            print(skipped)
        minimal_df = minimal_df.dropna()
        
        k_v = zip(
            minimal_df[key],
            minimal_df[value]
        )
        if minimal_df[key].duplicated().any():
            if coerce:
                warn(f'Reducing to a simple dictionary even though some keys'
                     f'point to more than one value (reason: coerce=True)')
            else:
                warn(f'Cannot reduce to simple dictionary: '
                     f'some keys point to more than one value; '
                     f'returned dict wil have sets of values instead; '
                     f'use coerce=True to override')
                d = {}
                for k, v in k_v:
                    if k not in d:
                        d[k] = set()
                    d[k].add(v)
                return d
            
        return dict(k_v)


class MetaSRA(MyDataFrame):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._lists_as_type('mapped_ontology_terms')
        self._lists_as_type('real-value_properties')
        self.accessions_map = SRAAccessionsMap()
    
    @classmethod
    def from_json(cls, path):
        meta_sra = read_json(path, orient='values')
        meta_sra = meta_sra.T
        meta_sra.columns = meta_sra.columns.str.replace(' ', '_')
        return MetaSRA(meta_sra)
    
    ontologies_data = {
        # uberon
        'anatomy': ('data/uberon/uberon.owl', identity),
        # disease_ontology
        'disease': ('data/disaese_ontology/HumanDO.owl', identity),
        # efo
        'experimental_factor': ('data/efo/efo.owl', identity),
        # cvcl
        'cellosaurus': ('data/cvcl/cellosaurus.obo', lambda i: i.replace(':', '_')),
        # cell ontology
        'cell': ('data/cl/cl.owl', identity)
    }
    
    _metadata = ['ontologies_data', 'accessions_map']
    
    def _lists_as_type(self, column, type=tuple):
        self[column] = self[column].apply(type)

    def expand_ontologies(self, to_strings=False):
        self = copy(self)
        for name, (ontology_path, id_mapper) in self.ontologies_data.items():
            ontology = Ontology(ontology_path)
            converter = (lambda entry: entry.name) if to_strings else identity
            self[name] = Series(
                tuple({
                    converter(ontology[id_mapper(ontology_id)])
                    for ontology_id in row.mapped_ontology_terms
                    if id_mapper(ontology_id) in ontology
                })
                for row in self.itertuples()
            ).values
        return MetaSRA(self)
    
    @cached_property
    def expanded(self):
        return self.expand_ontologies(to_strings=False)
    
    def summarize_ontologies(self, to_strings=False):
        summary = {}
        converter = (lambda entry: entry.name) if to_strings else identity
        for name in self.ontologies_data:
            counter = Counter()
            for terms_list in self.expanded[name]:
                counter.update([converter(term) for term in terms_list])
            counts = Series(dict(counter.most_common()))
            counts.name = 'count'
            summary[name] = DataFrame(counts)
        return HorizontalNamespace(summary)
    
    def limit_by_studies(self, studies):
        # sample = SRS
        # study = SRP
        srs_to_srp = self.accessions_map.sample_to_study(coerce=True)
        studies = set(studies)
        if 'study' not in self.columns:
            self['study'] = self.index.map(srs_to_srp)
        return MetaSRA(self[self.study.isin(studies)])
    
    @property
    def diseased(self):
        return self[self.mapped_ontology_terms.apply(
            lambda terms: any(
                term.startswith('DOID')
                for term in terms
            )
        )]
    
    def run_to_disease(self, run):
        try:
            sample = self.accessions_map.accession_to_sample()[run]
            return self.expanded[sample].disease_ontology
        except KeyError:
            return nan
    
    def get_ontologies(self):
        ontologies = set()
        self.mapped_ontology_terms.apply(
            lambda terms: ontologies.update([
                term.split(':')[0] for term in terms
            ])
        )
        return ontologies
    
    def get_diseases(self):
        diseases = set()
        self.mapped_ontology_terms.apply(
            lambda terms: diseases.update([
                term
                for term in terms
                if term.startswith('DOID:')
            ])
        )
        assert all(disease.startswith('DOID:') for disease in diseases)
        return diseases
    
    def limit_by_doid(self, diseases):
        diseases = set(diseases)
        return self[
            self.mapped_ontology_terms.apply(
                lambda terms: any(
                    term in diseases
                    for term in terms
                )
            )
        ]
