from pandas import DataFrame

from data_frames import MyDataFrame
from . import ExpressionManager


def get_subtype_by_sample(expression: ExpressionManager, subtype_by_participant: DataFrame, subtype_column):

    samples_by_participant = expression.samples_by_participant()
    participants = set(samples_by_participant)

    subtypes = MyDataFrame(subtype_by_participant[
       subtype_by_participant.participant.isin(participants)
    ])
    subtypes.normalize_columns()

    barcode_participant_sample = DataFrame([
        {'barcode': barcode, 'participant': participant, 'sample': barcode.barcode}
        for participant, samples in samples_by_participant.items()
        for barcode in samples
    ])

    merged = barcode_participant_sample.merge(subtypes)
    return merged[['sample', subtype_column]]


def group_by_subtype(subtype_sample_df, subtype_column) -> dict:
    """Caveat: a single participant can have multiple different tumour sub-types."""
    return subtype_sample_df.groupby(subtype_column)['sample'].apply(lambda samples: samples.tolist()).to_dict()

