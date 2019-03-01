from copy import copy
from io import StringIO
from typing import Mapping, Any

from pandas import read_table


class Barcode:
    pass


class TCGABarcode(Barcode):

    # Aliquot barcode is the highest resolution identifier we can get in TCGA, e.g.:
    # TCGA-02-0001-01B-02D-0182-06

    # copy-paste from https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
    annotations_table = """
    Label 	Identifier for 	Value 	Value Description 	Possible Values
    Analyte 	Molecular type of analyte for analysis 	D 	The analyte is a DNA sample 	See Code Tables Report
    Plate 	Order of plate in a sequence of 96-well plates 	182 	The 182nd plate 	4-digit alphanumeric value
    Portion 	Order of portion in a sequence of 100 - 120 mg sample portions 	1 	The first portion of the sample 	01-99
    Vial 	Order of sample in a sequence of samples 	C 	The third vial 	A to Z
    Project 	Project name 	TCGA 	TCGA project 	TCGA
    Sample 	Sample type 	1 	A solid tumor 	Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29. See Code Tables Report for a complete list of sample codes
    Center 	Sequencing or characterization center that will receive the aliquot for analysis 	1 	The Broad Institute GCC 	See Code Tables Report
    Participant 	Study participant 	1 	The first participant from MD Anderson for GBM study 	Any alpha-numeric value
    TSS 	Tissue source site 	2 	GBM (brain tumor) sample from MD Anderson 	See Code Tables Report
    """

    metadata = read_table(
        StringIO(annotations_table),
        delimiter=' *\t *',
        engine='python'
    )

    project: str
    tss: str
    participant: str
    sample: str
    vial: str
    portion: str
    analyte: str
    plate: str
    center: str

    dashed_parts_names = [
        'project',
        'tss',
        'participant',
        'sample_and_vial',
        'portion_and_analyte',
        'plate',
        'center'
    ]

    joined_alphanumeric_parts = {
        'sample_and_vial': [
            'sample',
            'vial',
        ],
        'portion_and_analyte': [
            'portion',
            'analyte'
        ]
    }

    def __init__(self, barcode):

        parsed_parts = self.decompose_barcode(barcode)

        for part_name, part_value in parsed_parts.items():
            setattr(self, part_name, part_value)

        self.barcode = barcode
        self.parsed_parts = parsed_parts

    def annotate(self):
        annotations = copy(self.metadata)
        annotations.rename(columns={'Value Description': 'Description of example value'}, inplace=True)
        annotations['Value'] = annotations['Label'].map(
            lambda label: getattr(self, label.lower())
        )
        return annotations

    def decompose_barcode(self, barcode: str) -> Mapping[str, Any]:
        dashed_parts = barcode.split('-')
        parsed_parts = dict(zip(self.dashed_parts_names, dashed_parts))
        for joined_name, (numeric_name, alpha_name) in self.joined_alphanumeric_parts.items():
            joined = parsed_parts.pop(joined_name, '')

            numeric_part = joined[:-1]
            alpha_part = joined[-1:]

            parsed_parts[numeric_name] = numeric_part
            parsed_parts[alpha_name] = alpha_part
        return parsed_parts

    def up_to(self, last_part: str) -> str:
        barcode_prefix = ''
        for part_name in self.dashed_parts_names:
            if part_name in self.joined_alphanumeric_parts:
                numeric, alpha = self.joined_alphanumeric_parts[part_name]
                if numeric != last_part:
                    barcode_prefix += '-' + self.parsed_parts[numeric]
                if alpha != last_part:
                    barcode_prefix += self.parsed_parts[alpha]
            else:
                barcode_prefix += '-' + self.parsed_parts[part_name]
            if part_name == last_part:
                break
        # remove the '-' at the beginning
        return barcode_prefix[1:]

    # From the metadata table:
    #     Tumor types range from 01-09, normal types from 10-19 and control samples from 20-29.
    #     See Code Tables Report for a complete list of sample codes.
    sample_type_ranges = {
        range(1, 10): 'tumor',
        range(11, 20): 'normal',
        range(21, 30): 'control',
    }

    @property
    def sample_type(self):
        for value_range, value_type in self.sample_type_ranges.items():
            if int(self.sample) in value_range:
                return value_type

    @classmethod
    def get_sample_type(cls, barcode: str) -> str:
        return TCGABarcode(barcode).sample_type

    @property
    def up_to_sample_type(self):
        parts = self.barcode.split('-')
        return TCGABarcode('-'.join(parts[:4])[:-1])

    def __hash__(self):
        return hash(self.barcode)

    def __repr__(self):
        return f'<TCGA Barcode {self.barcode}>'
    
    def __eq__(self, other):
        return self.barcode == other.barcode
