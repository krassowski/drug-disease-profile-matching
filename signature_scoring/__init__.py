from pandas import Series

from .models import SignaturesGrouping
from .processor import SignatureProcessor
from .processor.fold_change import FoldChangeSignatureProcessor


def score_signatures(
    scoring_func, disease_signature, signatures=None, limit=500, gene_subset=None,
    fold_changes=False, scale=False, progress=False, processes=None,
    gene_selection=Series.nlargest, warning_manager=None
):
    processor_type = SignatureProcessor

    if fold_changes:
        print('Make sure that you use disease_signature with fold changes as well')
        processor_type = FoldChangeSignatureProcessor

    processor = processor_type(
        signatures if isinstance(signatures, SignaturesGrouping) else scoring_func.collection(signatures),
        warning_manager,
        progress,
        processes
    )
    return processor.score_signatures(
        scoring_func, disease_signature, limit, gene_subset, scale, gene_selection
    )

