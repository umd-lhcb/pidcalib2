import pytest

from pidcalib2 import binning


def test_p_binning():
    with pytest.raises(KeyError):
        binning.p_binning("graviton")

    assert binning.p_binning("pi") == [
        3000,
        9300,
        15600,
        19000.0,
        24400.0,
        29800.0,
        35200.0,
        40600.0,
        46000.0,
        51400.0,
        56800.0,
        62200.0,
        67600.0,
        73000.0,
        78400.0,
        83800.0,
        89200.0,
        94600.0,
        100000.0,
    ]
