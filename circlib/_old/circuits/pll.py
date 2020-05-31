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
    Type-I digital phase detector circuit.
    """

    def __init__(self, vdd=3.3, rval=1e2, cval=1e-9):
        self.vdd = vdd
        self.rval = rval
        self.cval = cval

    @subcircuit
    def cad(self) -> Subcircuit:
        """
        """
        raise NotImplementedError("Cad not yet implemented.")

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


class Type2PhaseDetector(BaseCircuit):
    """
    Type 2 phase detector circuit.  This only works with digital
    signals.
    """

    def __init__(self, vdd=3.3):
        """
        """
        self.vdd = vdd
        self.pins = None

    def _assign_pins(self):
        self.pins = {
            "fref": self.adc["in"][0],
            "fvco": self.adc["in"][1],
            "vout": self.limit["out"],
            "vcc": self.isourcet["op"],
            "gnd": self.c["n"],
        }

    def __getitem__(self, key):
        return self.pins[key]

    def __setitem__(self, key, value):
        self.pins[key] = value

    @subcircuit
    def cad(self) -> Subcircuit:
        """
        """
        raise NotImplementedError("Cad not yet implemented.")

    @subcircuit
    def spice(self, current_gain=1e-3):
        """
        :returns: [fref, fvco, vout, vcc, gnd]
        """
        self.vp = Part("pyspice", "V", value=self.vdd)

        dff_model = XspiceModel(
            "dff",
            "d_dff",
            clk_delay=1e-12,
            set_delay=1e-12,
            reset_delay=1e-12,
            ic=0,
            data_load=1e-12,
            clk_load=1e-12,
            set_load=1e-12,
            reset_load=1e-12,
            rise_delay=1e-12,
            fall_delay=1e-12,
        )
        inv_model = XspiceModel(
            "inverter",
            "d_inverter",
            rise_delay=2e-9,
            fall_delay=2e-9,
            input_load=1e-12,
        )

        self.dt = Part(
            "pyspice",
            "A",
            io=["d", "clk", "null", "rst", "q", "null"],
            model=dff_model,
        )
        self.db = Part(
            "pyspice",
            "A",
            io=["d", "clk", "null", "rst", "q", "null"],
            model=dff_model,
        )
        self.inv1 = Part("pyspice", "A", io=["in", "out"], model=inv_model)
        self.inv2 = Part("pyspice", "A", io=["in", "out"], model=inv_model)
        self.and_gate = Part(
            "pyspice",
            "A",
            io=["in[]", "out"],
            model=XspiceModel(
                "and",
                "d_and",
                rise_delay=1e-12,
                fall_delay=1e-12,
                input_load=1e-12,
            ),
        )
        self.isourcet = Part("pyspice", "G", current_gain=current_gain)
        self.isourceb = Part("pyspice", "G", current_gain=current_gain)
        self.limit = Part(
            "pyspice",
            "A",
            io=["in", "out"],
            model=XspiceModel(
                "limiter",
                "limit",
                in_offset=0,
                gain=1,
                out_lower_limit=0,
                out_upper_limit=1,
            ),
        )
        self.c = Part("pyspice", "C", value=1e-9)
        self.r = Part("pyspice", "R", value=1e6)

        self.adc = Part(
            "pyspice",
            "A",
            io=["in[]", "out[]"],
            model=XspiceModel(
                "adc",
                "adc_bridge",
                in_low=self.vdd / 3,
                in_high=2 * self.vdd / 3,
                rise_delay=1e-12,
                fall_delay=1e-12,
            ),
        )
        self.dac = Part(
            "pyspice",
            "A",
            io=["in[]", "out[]"],
            model=XspiceModel(
                "dac",
                "dac_bridge",
                out_low=0,
                out_high=self.vdd,
                out_undef=self.vdd / 2,
                input_load=1e-12,
                t_rise=1e-12,
                t_fall=1e-12,
            ),
        )

        self._connect_components()
        self._assign_pins()

    def _connect_components(self):
        """
        """
        # ground net
        self.isourcet["in"] += (
            self.isourceb["on"],
            self.isourceb["in"],
            self.c["n"],
            self.vp["n"],
            self.r["n"],
        )

        self.dt["rst"] += self.db["rst"], self.inv2["out"]
        self.inv1["out"] += self.inv2["in"]
        self.inv1["in"] += self.and_gate["out"]
        self.dt["q"] += self.and_gate["in"][0], self.dac["in"][0]
        self.db["q"] += self.and_gate["in"][1], self.dac["in"][1]
        self.dac["out"][0] += self.isourcet["ip"]
        self.dac["out"][1] += self.isourceb["ip"]
        self.isourcet["on"] += (
            self.isourceb["op"],
            self.limit["in"],
            self.c["p"],
            self.r["p"],
        )
        self.adc["out"][0] += self.dt["clk"]
        self.adc["out"][1] += self.db["clk"]
        self.vp["p"] += self.adc["in"][2]
        self.adc["out"][2] += self.dt["d"], self.db["d"]
