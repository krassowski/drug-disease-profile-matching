from data_sources.tcga.barcode import TCGABarcode


def test_barcode_parts():

    barcode = TCGABarcode('TCGA-02-0001-01B-02D-0182-06')

    assert barcode.project == 'TCGA'

    assert barcode.tss == '02'

    assert barcode.participant == '0001'

    assert barcode.sample == '01'
    assert barcode.vial == 'B'

    assert barcode.portion == '02'
    assert barcode.analyte == 'D'

    assert barcode.plate == '0182'

    assert barcode.center == '06'

    assert barcode.sample_type == 'tumor'


def test_incomplete_barcode():

    barcode = TCGABarcode('TCGA-OR-A5J1')

    assert barcode.project == 'TCGA'

    assert barcode.tss == 'OR'

    assert barcode.participant == 'A5J1'
