import re
from glob import glob
from copy import copy

from pandas import read_csv, read_table
from tempfile import TemporaryDirectory, NamedTemporaryFile
from os import makedirs
from pathlib import Path
from warnings import warn

from config import DATA_DIR
from models import ExpressionProfile


class GSEA:

    def __init__(self):
        self.msigdb = MolecularSignaturesDatabase()

    def prepare_files(self, expression_data: ExpressionProfile):
        with NamedTemporaryFile('w', delete=False, suffix='.cls') as f:
            classes = expression_data.classes
            f.write(f'{len(classes)} {len(set(classes))} 1\n')
            classes_set = []
            for class_ in classes:
                if class_ not in classes_set:
                    classes_set.append(class_)
            f.write(f'# {" ".join(classes_set)}\n')
            f.write(' '.join(classes))
            classes = f.name

        with NamedTemporaryFile('w', delete=False, suffix='.txt') as f:
            expression_data = copy(expression_data)
            columns = expression_data.columns
            expression_data['DESCRIPTION'] = 'na'
            expression_data = expression_data[['DESCRIPTION', *columns]]
            if type(expression_data.index[0]) is bytes:
                expression_data.index = [b.decode('utf-8') for b in expression_data.index]
            expression_data.index.name = 'gene'
            expression_data.to_csv(f, sep='\t')
            data = f.name

        return data, classes

    def forward_streams(self, stdout, stderr):
        def display_one_at_a_time(text):
            from IPython.core.display import clear_output
            clear_output(wait=True)
            from IPython.core.display import display
            display(text)
        streams = {
            #stdout: display_one_at_a_time,
            stderr: print  # warn is silenced?
        }
        while True:
            active_streams = {}
            for stream, handler in streams.items():

                line = stream.readline()
                if line:
                    line = line.decode('utf-8')
                    handler(line.rstrip())
                    active_streams[stream] = handler

            streams = active_streams
            if not streams:
                return


class GeneSet:

    def __init__(self, name, genes):
        self.name = name
        self.genes = set(genes)

    @classmethod
    def from_gmt_line(cls, line):
        name, url, *ids = line.split()
        return cls(name, ids)


class MolecularSignaturesDatabase:
    def __init__(self, version='6.2'):
        self.path = Path(DATA_DIR + '/msigdb')
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

    def load(self, gene_sets, id_type):
        path = self.resolve(gene_sets=gene_sets, id_type=id_type)

        with open(path) as f:
            return {
                GeneSet.from_gmt_line(line)
                for line in f
            }


class GSEADesktop(GSEA):
    path = Path('java')

    def __init__(self, memory_size=5000):
        super().__init__()
        self.memory_size = memory_size

    def core_command(self):
        pwd = Path(__file__).absolute().parent.parent
        gsea_path = pwd / 'thirdparty' / 'gsea-3.0.jar'
        assert gsea_path.exists()
        return str(self.path) + f' -cp {str(gsea_path)}'

    def load_results(self, outdir, name, expression_data):
        timestamps = []
        for path in glob(f'{outdir}/{name}.Gsea.*'):
            match = re.match(f'{name}\.Gsea\.(\d*)', Path(path).name)
            timestamp = match.group(1)
            timestamps.append(timestamp)

        timestamps = sorted(timestamps, reverse=True)
        newest = timestamps[0]

        results = {}
        for class_type in set(expression_data.classes):
            result = read_table(f'{outdir}/{name}.Gsea.{newest}/gsea_report_for_{class_type}_{newest}.xls')
            assert all(result['Unnamed: 11'].isnull())
            result.drop('Unnamed: 11', axis='columns', inplace=True)
            result.columns = [column.lower().replace(' ', '_') for column in result.columns]
            result.drop('gs<br>_follow_link_to_msigdb', axis='columns', inplace=True)
            result.set_index('name', inplace=True)
            results[class_type] = result

        return results

    def run(
        self, expression_data, gene_sets, id_type='symbols', collapse=True, mode='Max_probe',
        permutations=1000, detailed_sets_analysis=False, permutation_type='phenotype', outdir=None,
        plot_top_x=0, name='my_analysis', metric='Signal2Noise', verbose=False, use_existing=False,
        normalization='meandiv',
        **kwargs
    ):
        """
        Memory in MB

        id_type:
            symbols
            entrez

        permutation_type:
            Gene_set

        metrics:
            Diff_of_Classes
            log2_Ratio_of_Classes
            Ratio_of_Classes
            tTest
            Signal2Noise
        """
        if use_existing:
            try:
                results = self.load_results(outdir, name, expression_data)
                print(
                    'Re-using pre-calculated results for GSEA'
                    ' (to prevent this, set use_existing=False or clear the temporary GSEA folder)'
                )
                return results
            except IndexError:
                pass

        command = self.core_command() + f' -Xmx{self.memory_size}m xtools.gsea.Gsea'

        for key, value in kwargs.items():
            command += f' -{key} {value}'

        if not gene_sets.endswith('.gmt'):
            gene_sets = self.msigdb.resolve(gene_sets, id_type)

        data_path, classes_path = self.prepare_files(expression_data)

        with TemporaryDirectory() as temp_dir:
            if not outdir:
                outdir = Path(temp_dir)
            else:
                outdir = Path(outdir)
                if not outdir.is_absolute():
                    pwd = Path().absolute().parent
                    outdir = pwd / str(outdir)
                makedirs(outdir, exist_ok=True)

            command += (
                f' -res {data_path}'
                f' -cls {classes_path}'
                f' -gmx {gene_sets}'
                f" -collapse false -mode {mode} -norm {normalization}"
                f" -nperm {permutations} -permute {permutation_type} -rnd_type no_balance"
                f" -scoring_scheme weighted -rpt_label {name}"
                f" -metric {metric} -sort real -order descending"
                " -create_gcts false -create_svgs false"
                # " -include_only_symbols true"
                f" -make_sets {str(detailed_sets_analysis).lower()} -median false -num 100"
                f" -plot_top_x {plot_top_x} -rnd_seed timestamp"
                " -save_rnd_lists false"
                " -set_max 500 -set_min 15"
                " -zip_report false"
                f' -out {outdir}'
                " -gui false"
            )
            if verbose:
                print(command)
            arguments = command.split(' ')
            from subprocess import PIPE, Popen
            process = Popen(arguments, stdout=PIPE, stderr=PIPE)
            self.forward_streams(process.stdout, process.stderr)

            results = self.load_results(outdir, name, expression_data)

            return results


class GSEApy(GSEA):
    path = Path('~/.pyenv/versions/drug_discovery/bin/gseapy')

    def __init__(self):
        super().__init__()
        import gseapy
        self.libraries = gseapy.get_library_name()

    def run(self, expression_data, program='gsea', id_type='symbols', gene_sets='KEGG_2016', outdir=None,
            permutation_type='phenotype', processes=8, verbose=True, plot=False, **kwargs):

        if not gene_sets.endswith('.gmt'):
            if gene_sets in self.libraries:
                assert id_type == 'symbols'
            else:
                gene_sets = self.msigdb.resolve(gene_sets, id_type)

        with TemporaryDirectory() as temp_dir:
            if not outdir:
                outdir = temp_dir
            outdir = Path(outdir)

            data_path, classes_path = self.prepare_files(expression_data)

            kwargs['data'] = data_path
            kwargs['cls'] = classes_path
            kwargs['permu-type'] = permutation_type
            kwargs['threads'] = processes
            kwargs['gmt'] = gene_sets
            kwargs['outdir'] = str(outdir)

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
            print(command)
            arguments = command.split(' ')
            from subprocess import PIPE, Popen
            process = Popen(arguments, stdout=PIPE, stderr=PIPE)
            self.forward_streams(process.stdout, process.stderr)
            result = read_csv(outdir / 'gseapy.gsea.phenotype.report.csv')
            return result
