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


class InvertingAmplifier(BaseCircuit):
    """
    Inverting amplifier circuit.
    """

    def __init__(self, r1_val, r2_val, amp_name="LF411C"):
        self.r1_val = r1_val
        self.r2_val = r2_val
        self.amp_name = amp_name
        lib_search_paths[SPICE].append("/home/matt/src/spicelib")

    @subcircuit
    def cad(self) -> Subcircuit:
        raise NotImplementedError()

    @subcircuit
    def spice(self) -> Subcircuit:
        """
        :returns: [in, out, vp, vm, gnd]
        """
        self.amp = Part(self.amp_name, self.amp_name)
        self.r1 = Part("pyspice", "R", value=self.r1_val)
        self.r2 = Part("pyspice", "R", value=self.r2_val)
        self._connect_components()
        return Subcircuit(
            pins=[
                self.r1["p"],
                self.r2["p"],
                self.amp["3"],
                self.amp["4"],
                self.amp["1"],
            ]
        )

    def _connect_components(self):
        # TODO amplifier should abstract interface to pins
        self.amp["2"] += self.r1["n"], self.r2["n"]
        self.amp["5"] += self.r2["p"]
