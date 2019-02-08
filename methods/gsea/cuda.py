from warnings import warn

from pandas import DataFrame

from config import third_party_dir
from helpers.streams import handle_streams
from signature_scoring.models.with_controls import ExpressionWithControls

from .base import GSEA
from methods.gsea.exceptions import GSEANoResults


class cudaGSEA(GSEA):
    path = third_party_dir / 'cudaGSEA/cudaGSEA/src/cudaGSEA'

    def __init__(self, fdr='full', use_cpu=False, **kwargs):
        super().__init__(**kwargs)
        assert self.path.exists()
        self.fdr = fdr
        self.log_tail = []
        self.use_cpu = use_cpu

    def read_from_process(self, process, handlers=None, verbose=False):

        results = []

        def save_result_or_display(text):
            if text.startswith('RESULT: '):
                statistics, gene_set = text[8:].split('\t')
                statistics = statistics.split(' ')
                result = {
                    statistics[i].rstrip(':'): float(statistics[i + 1])
                    for i in range(0, len(statistics), 2)
                }
                assert gene_set.startswith('(') and gene_set.endswith(')')
                result['gene_set'] = gene_set[1:-1]
                results.append(result)
            else:
                self.log_tail.append(text)
                if len(self.log_tail) > 4:
                    self.log_tail = self.log_tail[1:]
                if verbose:
                    self.display(text)

        handlers = {
            'out': save_result_or_display,
            'err': warn
        }

        handle_streams(process, handlers)

        return DataFrame(results)

    metrics_using_variance = {
        'onepass_t_test',
        'onepass_signal2noise',
        'twopass_signal2noise',
        'twopass_t_test',
        'stable_signal2noise',
        'stable_t_test',
        'overkill_signal2noise',
        'overkill_t_test'
    }

    def run(
        self, expression_data: ExpressionWithControls, gene_sets: str,
        metric='twopass_signal2noise', id_type='symbols',
        permutations=1000, permutation_type='phenotype',
        verbose=False, delete=True, **kwargs
    ):
        super().run(expression_data, gene_sets, metric=metric)

        assert permutation_type == 'phenotype'

        if not gene_sets.endswith('.gmt'):
            gene_sets = self.msigdb.resolve(gene_sets, id_type)

        data_path, classes_path = self.prepare_files(expression_data)

        command = (
            f'{self.path}'
            f' -res {data_path}'
            f' -cls {classes_path}'
            f' -gmx {gene_sets}'
            f' -nperm {permutations}'
            f' -metric {metric}'
            ' -order descending'
            f' -fdr {self.fdr}'
            + (' -cpu' if self.use_cpu else '')
        )
        results = self.run_in_subprocess(command, verbose=verbose)

        if delete:
            self.clean_up(data_path, classes_path)

        cuda_pitfalls_warning = (
            'There might be an error with CUDA (especially after waking the computer up from sleep); '
            'if you run your graphics server (X/Wayland) on a separate card, '
            'you can try using: sudo rmmod nvidia_uvm; sudo modprobe nvidia_uvm'
        )

        if results.empty:
            last_line = self.log_tail[-1]
            if 'error' in last_line.lower():
                if not verbose:
                    print('Command:', command)
                if 'Error: duplicate symbol' in last_line:
                    print(
                        'Duplicate symbol error may be caused by duplicated genes in index, or '
                        'by having an empty column in the data. Please inspect this data frame: '
                    )
                    print(expression_data.joined)
                if 'CUDA error' in last_line:
                    warn(cuda_pitfalls_warning)
            print('Tail of the logs:')
            print(self.log_tail)
            raise GSEANoResults('No results returned.')

        results = results.set_index('gene_set').rename({
            'ES': 'es',
            'FWER': 'fwer_p-val',
            'NES': 'nes',
            'NP': 'nom_p-val',
            'FDR': 'fdr_q-val'
        }, axis=1)

        return {
            expression_data.case_name: results[results.nes > 0],
            expression_data.control_name: results[results.nes < 0]
        }
