#!/usr/bin/env python

from skidl import Part, subcircuit
from skidl.pyspice import gnd, node, generate_netlist, lib_search_paths, SPICE
from libcircuit.cad import (
    eseries_val,
    resistor_footprint,
    capacitor_footprint,
)
from libcircuit.spice import capacitor_equiv, resistor_equiv
import matplotlib.pyplot as plt


class CapacitanceMultiplier:
    """
    Capacitance multiplier circuit.
    """

    def __init__(self, iloadmax, ripplefreq, vin, vout, gnd):
        """
        Capacitance multiplier subcircuit. @iloadmax is the maximum load
        current of all downstream devices. @ripplefreq is the ripple
        frequency of the input voltage. When placing this at the output of
        a voltage converter, use the converter's switching frequency. @vin
        and @vout are the input and output nets, respectively.
        """
        self.rval = eseries_val(100 / iloadmax)
        self.cval = eseries_val(10 / (self.rval * ripplefreq))
        self.vin = vin
        self.vout = vout
        self.gnd = gnd

    @subcircuit
    def cad(self):
        """
        Generate a physical circuit for use in a CAD netlist.
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

    @subcircuit
    def spice(self, use_parasitics: bool = True):
        """
        Generate a spice netlist for simulation.
        """
        if use_parasitics:
            self.r = resistor_equiv(value=self.rval)
            self.c = capacitor_equiv(value=self.cval)
            self.rbase = resistor_equiv(value=100)
        else:
            self.r = Part("pyspice", "R", value=self.rval)
            self.c = Part("pyspice", "C", value=self.cval)
            self.rbase = Part("pyspice", "R", value=100)

        self.d = Part("pyspice", "D", model="1N4001")
        lib_search_paths[SPICE].append("/home/matt/src/spicelib")
        self.nmos = Part("FQD13N06", "FQD13N06")
        # self.npn = Part("pyspice", "Q", model="2N2222A")
        self._connect_components()

    def _connect_components(self):
        # self.vin += self.r[1], self.npn["C"]
        self.vin += self.r[1], self.nmos[2]
        # self.r[2] += self.c[1], self.rbase[2], self.d[2]
        self.r[2] += self.c[1], self.nmos[1], self.d[2]
        # self.rbase[1] += self.npn["B"]
        # self.rbase[1] += self.nmos[2]
        # self.npn["E"] += self.vout, self.d[1]
        self.nmos[3] += self.vout, self.d[1]
        self.c[2] += self.gnd


if __name__ == "__main__":
    # def simulate():
    vin = Part("pyspice", "SINEV", offset=5, amplitude=50e-3)
    rload = Part("pyspice", "R", value=10e3)
    rload[1]
    vin[2] += gnd, rload[2]
    cap_mul = CapacitanceMultiplier(
        iloadmax=1, ripplefreq=500e3, vin=vin[1], vout=rload[1], gnd=gnd
    )
    print(cap_mul.rval)
    print(cap_mul.cval)
    cap_mul.spice()
    circ = generate_netlist(libs="/home/matt/src/spicelib")
    sim = circ.simulator()
    waveforms = sim.transient(step_time=100e-6, end_time=100e-3)
    time = waveforms.time
    vin = waveforms[node(vin[1])]
    vout = waveforms[node(rload[1])]
    fig = plt.figure()
    plt.plot(time, vin)
    plt.plot(time, vout)
    plt.show()
    # print("time\tvin\tvout")
    # for (t, vi, vo) in zip(
    #     time.as_ndarray(), vin.as_ndarray(), vout.as_ndarray()
    # ):
    #     print("{}\t{}\t{}".format(t, vi, vo))
