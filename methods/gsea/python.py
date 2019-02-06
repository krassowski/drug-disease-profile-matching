from pathlib import Path
from tempfile import TemporaryDirectory

from pandas import read_csv

from config import third_party_dir

from .base import GSEA


class GSEApy(GSEA):
    """
    To use gsea.py please install it with:
        >>> pip3 install gseapy
    and link its binary to the thirdparty directory
        >>> ln -s virtual_environment_path/bin/gseapy thirdparty/gseapy
    """

    path = Path(third_party_dir / 'gseapy')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if not self.path.exists():
            print(self.__doc__)

        import gseapy
        self.libraries = gseapy.get_library_name()

    metrics_using_variance = {
        'signal_to_noise',
        't_test'
    }

    def run(
        self, expression_data, gene_sets='KEGG_2016', metric='signal_to_noise', id_type='symbols',
        permutations=1000, permutation_type='phenotype', verbose=True,
        program='gsea', processes=8, plot=False, out_dir=None, **kwargs
    ):

        super().run(expression_data, gene_sets, metric=metric)

        if not gene_sets.endswith('.gmt'):
            if gene_sets in self.libraries:
                assert id_type == 'symbols'
            else:
                gene_sets = self.msigdb.resolve(gene_sets, id_type)

        with TemporaryDirectory() as temp_dir:
            out_dir = self.prepare_out_dir(out_dir, temp_dir)

            data_path, classes_path = self.prepare_files(expression_data)

            kwargs.update({
                'data': data_path,
                'cls': classes_path,
                'permu-type': permutation_type,
                'permu-num': permutations,
                'threads': processes,
                'gmt': gene_sets,
                'method': metric,
                'outdir': str(out_dir)
            })

            command = ' '.join([
                str(self.path.expanduser()),
                program,
                *[
                    f'--{name} {value}'
                    for name, value in kwargs.items()
                ]
            ])
            command += ' --verbose' if verbose else ''
            command += ' --no-plot' if not plot else ''

            self.run_in_subprocess(command, verbose=verbose)

            result = read_csv(out_dir / f'gseapy.gsea.{permutation_type}.report.csv')
            return result
