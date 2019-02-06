import re
from glob import glob
from pathlib import Path
from typing import Set

from config import DATA_DIR


class GeneSet:

    def __init__(self, name, genes):
        self.name = name
        self.genes = set(genes)

    @classmethod
    def from_gmt_line(cls, line):
        name, url, *ids = line.split()
        return cls(name, ids)


class GeneMatrixTransposed:

    def __init__(self, gene_sets):
        self.gene_sets = gene_sets

    @classmethod
    def from_gmt(cls, path):
        with open(path) as f:
            return cls({
                GeneSet.from_gmt_line(line)
                for line in f
            })

    def trim(self, min_genes, max_genes: int):
        return GeneMatrixTransposed({
            gene_set
            for gene_set in self.gene_sets
            if min_genes <= len(gene_set.genes) <= max_genes
        })

    def to_gmt(self, path):
        with open(path, mode='w') as f:
            for gene_set in self.gene_sets:
                f.write(gene_set.name + '\t' + '\t'.join(gene_set.genes) + '\n')

    def subset(self, genes: Set[str]):
        return GeneMatrixTransposed({
            GeneSet(name=gene_set.name, genes=gene_set.genes & genes)
            for gene_set in self.gene_sets
        })


class MolecularSignaturesDatabase:
    def __init__(self, version='6.2'):
        self.path = Path(DATA_DIR) / 'msigdb'
        self.version = version
        wildcard_path = str((self.path / f'*.v{self.version}.*.gmt').resolve())
        self.gene_sets = [
            self.parse_name(Path(p).name)
            for p in glob(wildcard_path)
        ]

    def parse_name(self, name):
        parsed = re.match(rf'(?P<name>.*?)\.v{self.version}\.(?P<id_type>(entrez|symbols)).gmt', name)
        return parsed.groupdict()

    def resolve(self, gene_sets, id_type):
        path = self.path / f'{gene_sets}.v6.2.{id_type}.gmt'
        if path.exists():
            return str(path)
        else:
            raise ValueError('Unknown library!')

    def load(self, gene_sets, id_type) -> GeneMatrixTransposed:
        path = self.resolve(gene_sets=gene_sets, id_type=id_type)

        return GeneMatrixTransposed.from_gmt(path)
