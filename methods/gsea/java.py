import re
from glob import glob
from pathlib import Path
from shutil import rmtree
from tempfile import TemporaryDirectory

from IPython.core.display import display
from pandas import read_table, DataFrame

from config import third_party_dir

from .base import GSEA
from .paths import tmp_dir


class GSEADesktop(GSEA):
    path = Path('java')

    def __init__(self, memory_size=5000, **kwargs):
        super().__init__(**kwargs)
        self.memory_size = memory_size
        self.gsea_path = third_party_dir / 'gsea-3.0.jar'
        assert self.gsea_path.exists()

    @property
    def core_command(self):
        return f'{self.path} -cp {self.gsea_path}'

    @staticmethod
    def load_results(out_dir, name, expression_data, delete=True):
        timestamps = []
        for path in glob(f'{out_dir}/{name}.Gsea.*'):
            match = re.match(rf'{name}\.Gsea\.(\d*)', Path(path).name)
            timestamp = match.group(1)
            timestamps.append(timestamp)

        timestamps = sorted(timestamps, reverse=True)
        newest = timestamps[0]

        results = {}
        for class_type in set(expression_data.classes):
            path = f'{out_dir}/{name}.Gsea.{newest}/gsea_report_for_{class_type}_{newest}.xls'
            result: DataFrame = read_table(path)
            assert all(result['Unnamed: 11'].isnull())
            result.drop('Unnamed: 11', axis='columns', inplace=True)
            result.columns = [column.lower().replace(' ', '_') for column in result.columns]
            result.drop('gs<br>_follow_link_to_msigdb', axis='columns', inplace=True)
            result.set_index('name', inplace=True)
            results[class_type] = result.dropna(subset=['nes', 'fdr_q-val'])

        if delete:
            rmtree(f'{out_dir}/{name}.Gsea.{newest}', ignore_errors=True)

        try:
            assert (results[expression_data.case_name].nes > 0).all()
            assert (results[expression_data.control_name].nes < 0).all()
        except AssertionError:
            display(results[expression_data.case_name])
            display(results[expression_data.control_name])
            raise

        return results

    metrics_using_variance = {'tTest', 'Signal2Noise'}

    def run(
        # Common GSEA interface
        self, expression_data, gene_sets, metric='Signal2Noise', id_type='symbols',
        permutations=1000, permutation_type='phenotype', verbose=False,
        # GSEADesktop-specific
        collapse=True, mode='Max_probe', detailed_sets_analysis=False, out_dir=None,
        plot_top_x=0, name='my_analysis', normalization='meandiv', min_genes=15, max_genes=500,
        # internal, won't be passed to GSEA
        use_existing=False, delete=True,
        # additional arguments to be passed
        **kwargs
    ):
        """Memory in MB"""
        super().run(expression_data, gene_sets, metric=metric)

        assert id_type in {'symbols', 'entrez'}
        assert permutation_type in {'Gene_set', 'phenotype'}
        assert metric in {'Diff_of_Classes', 'log2_Ratio_of_Classes', 'Ratio_of_Classes' 'tTest', 'Signal2Noise'}

        if use_existing:
            try:
                results = self.load_results(out_dir, name, expression_data)
                print(
                    'Re-using pre-calculated results for GSEA'
                    ' (to prevent this, set use_existing=False or clear the temporary GSEA folder)'
                )
                return results
            except IndexError:
                pass

        command = f'{self.core_command} -Xmx{self.memory_size}m xtools.gsea.Gsea'

        for key, value in kwargs.items():
            command += f' -{key} {value}'

        if not gene_sets.endswith('.gmt'):
            gene_sets = self.msigdb.resolve(gene_sets, id_type)

        data_path, classes_path = self.prepare_files(expression_data)

        with TemporaryDirectory(dir=tmp_dir) as temp_dir:
            out_dir = self.prepare_out_dir(out_dir, temp_dir)

            command += (
                f' -res {data_path}'
                f' -cls {classes_path}'
                f' -gmx {gene_sets}'
                f' -collapse false -mode {mode} -norm {normalization}'
                f' -nperm {permutations} -permute {permutation_type} -rnd_type no_balance'
                f' -scoring_scheme weighted -rpt_label {name}'
                f' -metric {metric}'
                ' -sort real -order descending'
                ' -create_gcts false -create_svgs false'
                # ' -include_only_symbols true'
                f' -make_sets {str(detailed_sets_analysis).lower()} -median false -num 100'
                f' -plot_top_x {plot_top_x} -rnd_seed timestamp'
                ' -save_rnd_lists false'
                f' -set_max {max_genes} -set_min {min_genes}'
                ' -zip_report false'
                f' -out {out_dir}'
                ' -gui false'
            )
            self.run_in_subprocess(command, verbose=verbose)

            results = self.load_results(out_dir, name, expression_data, delete=delete)

        if delete:
            self.clean_up(data_path, classes_path)

        return results
