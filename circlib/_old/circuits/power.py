#!/usr/bin/env python

from skidl import Part, subcircuit
from skidl.pyspice import lib_search_paths, SPICE
from libcircuit.cad import (
    eseries_val,
    resistor_footprint,
    capacitor_footprint,
)
from libcircuit.base import BaseCircuit
from libcircuit.skidl import Subcircuit
from libcircuit.spice import capacitor_equiv, resistor_equiv


class CapacitanceMultiplier(BaseCircuit):
    """
    Capacitance multiplier circuit.
    """

    def __init__(self, iloadmax, ripplefreq):
        """
        Capacitance multiplier subcircuit.

        :param iloadmax: maximum load current of all downstream
            devices
        :param ripplefreq: ripple frequency of the input voltage.
            When placing this at the output of a voltage converter,
            use the converter's switching frequency.
        """
        self.rval = eseries_val(10e3 / iloadmax)
        self.cval = eseries_val(10 / (self.rval * ripplefreq))

    @subcircuit
    def cad(self) -> Subcircuit:
        """
        Generate a physical circuit for use in a CAD netlist.

        :returns: Subcircuit instance with pins [vin, vout, gnd]
        """
        self.r = Part(
            "Device",
            "R",
            value="{}".format(self.rval),
            footprint=resistor_footprint(self.rval),
        )
        self.c = Part(
            "Device",
            "C",
            value="{}".format(self.cval),
            footprint=capacitor_footprint(self.cval),
        )
        self.rbase = Part(
            "device", "R", value="100", footprint=resistor_footprint(100)
        )
        self.d = Part("Diode", "1N4001")
        # TODO should use a more "intelligent" selection process
        self.npn = Part("Transistor_BJT", "2STN1550")
        self._connect_components()
        return Subcircuit(pins=[self.r[1], self.transistor[3], self.c[2]])

    @subcircuit
    def spice(self, use_parasitics: bool = True) -> Subcircuit:
        """
        Generate a spice netlist for simulation.

        :param use_parasitics: Use non-ideal, parasitic models for
            resistors and capacitors.
        """
        if use_parasitics:
            self.r = resistor_equiv(value=self.rval)
            self.c = capacitor_equiv(value=self.cval)
            # self.rbase = resistor_equiv(value=100)
        else:
            self.r = Part("pyspice", "R", value=self.rval)
            self.c = Part("pyspice", "C", value=self.cval)
            self.rbase = Part("pyspice", "R", value=100)

        self.d = Part("pyspice", "D", model="1N4001")
        lib_search_paths[SPICE].append("/home/matt/src/spicelib")
        # self.transistor = Part("FQD13N06", "FQD13N06")
        self.transistor = Part("pyspice", "Q", model="2N2222A")
        self._connect_components()
        return Subcircuit(pins=[self.r[1], self.transistor[3], self.c[2]])

    def _connect_components(self):
        # self.vin += self.r[1], self.npn["C"]
        self.r[1] += self.transistor[2]
        # self.r[2] += self.c[1], self.rbase[2], self.d[2]
        self.r[2] += self.c[1], self.transistor[1], self.d[2]
        # self.rbase[1] += self.npn["B"]
        # self.rbase[1] += self.transistor[2]
        # self.npn["E"] += self.vout, self.d[1]
        self.transistor[3] += self.d[1]
