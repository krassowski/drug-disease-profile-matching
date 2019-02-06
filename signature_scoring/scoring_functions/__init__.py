from types import FunctionType
from typing import Type
from dataclasses import dataclass

from ..models import (
    Profile,
    SignaturesGrouping, SignaturesCollection
)
from signature_scoring.models.with_controls import ExpressionWithControls, SubstancesCollectionWithControls


@dataclass
class ScoringFunction:
    func: FunctionType

    # whether the query and signatures dataset should include:
    #   - case-biological_control array (ExpressionWithControls), or
    #   - profiles derived from precomputed differential signatures (Profile)
    #   - case-plate_control (NotImplemented, requires level4 data)
    # Rationale:
    #   case-control pairs might be required if a scoring function
    #   calculates differential signatures using it's own methods
    #   (e.g. mroast). In all other cases, use of precomputed
    #   signatures is more effective (as we don't have to recompute
    #   the differential profile for each call of the scoring function)
    input: Type = Profile

    # is grouping of signatures required by the scoring function?
    #   None (default) - each signature should be processed and passed to the scoring function individually
    #   'by_replicate' - not implemented yet, requires downloading level4 data
    #   'by_substance'
    # Rationale:
    #   mroast needs to estimate variance to compute differential profiles;
    #   therefore it requires signatures to be grouped by replicates.
    #   as for this versions, replicates handling is not implemented
    #   (as I do not have enough space to download it), so a trick of using
    #   per-substance grouping (i.e. same substance, different concentrations)
    #   was used to validate correctness of the pipeline.
    grouping: str = None

    custom_multiprocessing: bool = False

    multiprocessing_exclude_first: bool = True

    # a hook to run janitorial tasks (garbage collection, display preparation)
    # before start of a larger batch of tasks; completely optional
    before_batch: FunctionType = lambda: None

    # if a function supports results caching then
    # its func should accept additional keyword argument:
    #   warn_about_cache: bool
    supports_cache: bool = None

    @property
    def collection(self) -> Type[SignaturesGrouping]:
        """Provides constructor which (when applied to SignaturesData)
        will create a structure yielding data of required input type
        and of required grouping.

        Not all combinations of input/grouping are implemented/allowed.
        """
        available_collections = {
            (ExpressionWithControls, 'by_substance'): SubstancesCollectionWithControls.from_signatures,
            (Profile, None): SignaturesCollection,
        }
        return available_collections[(self.input, self.grouping)]

    @property
    def is_applicable_to_control_signatures(self):
        return self.input != ExpressionWithControls

    def __call__(self, disease: input, compound: input, **kwargs):
        return self.func(disease, compound, **kwargs)


def scoring_function(func, **kwargs):
    proxy = ScoringFunction(func, **kwargs)
    proxy.__name__ = func.__name__
    proxy.original_function = func
    if hasattr(func, '__source__'):
        proxy.__source__ = func.__source__
    return proxy


class ScoringError(Exception):
    pass
