from pandas import DataFrame, concat

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

    # try to merge by sample
    merged_by_sample = barcode_participant_sample.merge(subtypes, left_on='sample', right_on='pan.samplesid')

    print(f'{len(merged_by_sample)} matched exactly on sample ID')

    # ignore participants that were merged above
    merged_participants = set(merged_by_sample.participant_x)
    remaining_subtypes = subtypes[~subtypes.participant.isin(merged_participants)]
    remaining_bps = barcode_participant_sample[~barcode_participant_sample.participant.isin(merged_participants)]

    results = [merged_by_sample]

    if len(remaining_bps):
        print('Not all samples were matched by sample, attempting to match by participant')
        merged_by_participant = remaining_bps.merge(remaining_subtypes, on='participant')
        print(f'{len(merged_by_participant)} matched by participant ID')
        results.append(merged_by_participant)

    merged = concat(results)

    print(f'{len(barcode_participant_sample) - len(merged)} not matched')

    return merged[['sample', subtype_column]]


def group_by_subtype(subtype_sample_df, subtype_column) -> dict:
    """Caveat: a single participant can have multiple different tumour sub-types."""
    return subtype_sample_df.groupby(subtype_column)['sample'].apply(lambda samples: samples.tolist()).to_dict()

