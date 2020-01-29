from skidl import Part, subcircuit
from skidl.pyspice import lib_search_paths, SPICE, A, XspiceModel
from libcircuit.cad import (
    eseries_val,
    resistor_footprint,
    capacitor_footprint,
)
from libcircuit.base import BaseCircuit
from libcircuit.skidl import Subcircuit
from libcircuit.spice import capacitor_equiv, resistor_equiv


class Type1PhaseDetectorDigital(BaseCircuit):
    """
    Type-I phase detector circuit.
    """

    def __init__(self, vdd=3.3, rval=1e2, cval=1e-9):
        self.vdd = vdd
        self.rval = rval
        self.cval = cval

    @subcircuit
    def cad(self) -> Subcircuit:
        """
        """
        raise ValueError("Cad not yet implemented.")

    @subcircuit
    def spice(self, use_parasitics: bool = True) -> Subcircuit:
        """
        Generate a spice netlist for simulation.

        :param use_parasitics: Use non-ideal, parasitic models for
            resistors and capacitors.

        :returns: [in1, in2, out, gnd]
        """
        self.low_threshold = self.vdd / 2 - 0.5
        self.high_threshold = self.vdd / 2 + 0.5

        if use_parasitics:
            self.r = resistor_equiv(value=self.rval)
            self.c = capacitor_equiv(value=self.cval)
        else:
            self.r = Part("pyspice", "R", value=self.rval)
            self.c = Part("pyspice", "C", value=self.cval)

        self.adc_model = XspiceModel(
            "adc",
            "adc_bridge",
            in_low=self.low_threshold,
            in_high=self.high_threshold,
            rise_delay=1e-12,
            fall_delay=1e-12,
        )
        self.adc = Part(
            "pyspice", "A", io=["anlg_in[]", "dig_out[]"], model=self.adc_model
        )
        self.dac = Part(
            "pyspice",
            "A",
            io=["dig_in[]", "anlg_out[]"],
            model=XspiceModel(
                "dac", "dac_bridge", out_low=0.0, out_high=self.vdd
            ),
        )
        self.xor = Part(
            "pyspice",
            "A",
            io=["in[]", "out"],
            model=XspiceModel(
                "xor",
                "d_xor",
                rise_delay=1e-12,
                fall_delay=1e-12,
                input_load=1e-12,
            ),
        )
        self._connect_components()
        return Subcircuit(
            pins=[
                self.adc["anlg_in"][0],
                self.adc["anlg_in"][1],
                self.c[1],
                self.c[2],
            ]
        )

    def _connect_components(self):
        """
        """
        self.adc["dig_out"][0] += self.xor["in"][0]
        self.adc["dig_out"][1] += self.xor["in"][1]
        self.xor["out"] += self.dac["dig_in"][0]
        self.dac["anlg_out"][0] += self.r[1]
        self.r[2] += self.c[1]
