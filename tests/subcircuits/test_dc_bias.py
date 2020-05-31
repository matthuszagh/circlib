import pytest
import numpy as np
from circlib.subcircuits.dc_bias import ResistiveDivider
from circlib.subcircuits.resistor import Resistor
from circlib.quantity import DCVoltage, DCCurrent
from circlib.eseries import Eseries


def test_rt_ubound1():
    div = ResistiveDivider(
        vp=DCVoltage(typval=10, minval=9.8, maxval=10.2),
        vn=DCVoltage(0),
        vout=DCVoltage(typval=3.3, minval=3.1, maxval=3.5),
        iout=DCCurrent(typval=1.5e-6, minval=1e-6, maxval=2e-6),
        resistor=Resistor(
            eseries=Eseries(series=24, tol=1e-2, limits=(1, 10e6))
        ),
    )
    assert np.isclose(div._rt_ubound(), 542857.1428571418)


def test_rt_ubound2():
    div = ResistiveDivider(
        vp=DCVoltage(typval=10, minval=9.8, maxval=10.2),
        vn=DCVoltage(0),
        vout=DCVoltage(typval=3.3, minval=3.1, maxval=3.5),
        iout=DCCurrent(typval=1.5e-3, minval=1e-3, maxval=2e-3),
        resistor=Resistor(
            eseries=Eseries(series=24, tol=1e-2, limits=(1, 10e6))
        ),
    )
    assert np.isclose(div._rt_ubound(), 542.8571428571418)


def test_set_top_bottom_singles1():
    div = ResistiveDivider(
        vp=DCVoltage(typval=10, minval=9.8, maxval=10.2),
        vn=DCVoltage(0),
        vout=DCVoltage(typval=3.3, minval=3.1, maxval=3.5),
        iout=DCCurrent(typval=1.5e-6, minval=1e-6, maxval=2e-6),
        resistor=Resistor(
            eseries=Eseries(series=24, tol=1e-2, limits=(1, 10e6))
        ),
    )
    assert np.isclose(div.r_top.resistance().typval, 47e3)
    assert np.isclose(div.r_bot.resistance().typval, 24e3)
    assert np.isclose(div.vout().typval, 3.309781690140844)
    assert np.isclose(div.vout().maxval, 3.4471537643069095)
    assert np.isclose(div.vout().minval, 3.174019707988208)


def test_set_top_bottom_singles2():
    div = ResistiveDivider(
        vp=DCVoltage(typval=10, minval=9.8, maxval=10.2),
        vn=DCVoltage(0),
        vout=DCVoltage(typval=3.3, minval=3.1, maxval=3.5),
        iout=DCCurrent(typval=1.5e-3, minval=1e-3, maxval=2e-3),
        resistor=Resistor(
            eseries=Eseries(series=24, tol=1e-2, limits=(1, 10e6))
        ),
    )
    assert np.isclose(div.r_top.resistance().typval, 47)
    assert np.isclose(div.r_bot.resistance().typval, 24)
    assert np.isclose(div.vout().typval, 3.309781690140844)
    assert np.isclose(div.vout().maxval, 3.4471537643069095)
    assert np.isclose(div.vout().minval, 3.174019707988208)
