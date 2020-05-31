import pytest
from circlib.quantity import Quantity, evaluate


def test_init():
    q = Quantity(8.1)
    assert q.minval == 8.1
    assert q.typval == 8.1
    assert q.maxval == 8.1

    q = Quantity(typval=5.5, minval=1)
    assert q.minval == 1
    assert q.typval == 5.5
    assert q.maxval == 5.5

    q = Quantity(minval=1, maxval=10.1)
    assert q.minval == 1
    assert q.typval == 11.1 / 2
    assert q.maxval == 10.1

    with pytest.raises(ValueError):
        q = Quantity()


def test_evaluate():
    """
    """
    a = Quantity(minval=8, maxval=10, typval=9)
    b = Quantity(minval=1, maxval=5, typval=4)
    res = evaluate((a - b) / a)
    assert res.minval == (a.minval - b.maxval) / a.minval
    assert res.maxval == (a.maxval - b.minval) / a.maxval
    assert res.typval == (a.typval - b.typval) / a.typval
