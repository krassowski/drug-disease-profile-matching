from functools import lru_cache
from typing import Tuple

import numpy as np
from pandas import DataFrame, Series

from data_sources.drug_connectivity_map import dcm

from .processor import progress_bar
from signature_scoring.profile import Signature
from .processor import SignatureProcessor


class FoldChangeSignature(Signature):
    def split(self, limit: int, nlargest=Series.nlargest) -> Tuple[Series, Series]:
        up_regulated = self[(self > 1) | (self > -1)]
        up_regulated = up_regulated[nlargest(up_regulated.abs(), limit).index]
        down_regulated = self[(self < 1) & (self < -1)]
        down_regulated = down_regulated[nlargest(-down_regulated.abs(), limit).index]
        return down_regulated, up_regulated


class FoldChangeSignatureProcessor(SignatureProcessor):

    signature = FoldChangeSignature

    def calculate_fold_change(self, signature, selected_genes, preserve_sign=False):

        controls = self.controls(selected_genes)
        control = controls[signature.name]

        problematic_genes = control[control == 0]
        if problematic_genes.any():
            print(problematic_genes)

        assert (control.index == signature.index).all()

        signature_fc = signature.divide(control)

        if preserve_sign:
            signs = signature.apply(np.sign)
            signature_fc = signature_fc.multiply(signs)

        return signature_fc

    def transform_signature(self, signature, selected_genes):
        return self.calculate_fold_change(signature, selected_genes)

    def controls(self, selected_genes):
        # sorted tuples to get hashable inputs for caching
        return self.get_controls_for_signatures(
            tuple(sorted(self.ids)),
            tuple(sorted(selected_genes))
        )

    @lru_cache()
    def get_controls_for_signatures(self, ids, genes_to_keep):
        controls_by_signature = {}
        for signature_id in progress_bar(ids):
            controls = dcm.get_controls(signature_id, exemplar_only=True)
            rows_to_keep = controls.index.isin(genes_to_keep)
            controls = controls[rows_to_keep]
            control = controls.mean(axis=1)
            controls_by_signature[signature_id] = control
        return DataFrame(controls_by_signature)
