from functools import lru_cache
from statistics import mean

from pandas import read_table
from pyfaidx import Fasta

from config import DATA_DIR


class ProteinSequences:
    """Implemented with Uniprot"""

    def __init__(self, file):

        def extract_id(header):
            return header.split('|')[1]

        self.fasta = Fasta(file, key_function=extract_id)

        id_mapping = read_table(
            DATA_DIR + '/uniprot/HUMAN_9606_idmapping.dat.gz',
            names=['uniprot_id', 'type', 'value']
        )
        gene_mappings = id_mapping[
            id_mapping.type.isin(['Gene_Name', 'Gene_Synonym'])
        ]
        self.gene_to_uniprot = gene_mappings[['value', 'uniprot_id']].set_index('value')

    @lru_cache()
    def average_length(self, gene_name=None):
        if gene_name:
            sequences = [
                self.fasta[uniprot_id]
                for uniprot_id in self.gene_to_uniprot.loc[gene_name].uniprot_id
                if uniprot_id in self.fasta
            ]
        else:
            sequences = self.fasta
        # TODO: warn if len > 2 and difference is significant
        return mean([len(sequence) for sequence in sequences])
