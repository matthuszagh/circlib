from skidl import subcircuit, Part
import numpy as np
from libcircuit.skidl import Subcircuit


def df_to_esr(c, df, df_freq):
    """
    Returns the equivalent series resistance (ESR) of a capacitor from
    a given dissipation factor, @df. @df_freq is the frequency at
    which DF was measured in hertz. @c is the capacitor's capacitance
    value in farads.
    """
    return df / (2 * np.pi * df_freq * c)


@subcircuit
def capacitor_equiv(
    value: float,
    rp: float = 500e6,
    rs: float = 50e-3,
    l: float = None,
    srf: float = 100e6,
) -> Subcircuit:
    """
    Non-ideal capacitor equivalent circuit for use in a SPICE
    simulation.

    Args:
        value: capacitance value (farads)
        rp: parallel (insulation) resistance (ohms)
        rs: ESR (ohms)
        l: series inductance (henrys)

    Returns:
        Subcircuit consisting of attachment pins.
    """

    if l is None:
        if srf is None:
            raise ValueError(
                "Must either specify inductance or self-resonant frequency "
                "for parasitic capacitor."
            )
        l = 1 / (((2 * np.pi * srf) ** 2) * value)

    # ============================= parts ============================
    rseries = Part("pyspice", "R", value=rs)
    rparallel = Part("pyspice", "R", value=rp)
    lseries = Part("pyspice", "L", value=l)
    cseries = Part("pyspice", "C", value=value)

    # ========================== connections =========================
    rseries[1] += rparallel[1]
    rseries[2] += lseries[1]
    lseries[2] += cseries[1]
    cseries[2] += rparallel[2]
    return Subcircuit(pins=[rseries[1], cseries[2]])


# TODO are these defaults reasonable? I'm having trouble determining
# resistor parasitics.
@subcircuit
def resistor_equiv(
    value: float, c: float = 0.4e-12, l: float = 2e-9
) -> Subcircuit:
    """
    Non-ideal resistor equivalent circuit for use in a SPICE
    simulation. @r is the resistor's resistance value (in
    ohms). @c is the parallel capacitance in ohms. @l is the series
    inductance in henrys.
    """
    # ============================= parts ============================
    rseries = Part("pyspice", "R", value=value)
    cparallel = Part("pyspice", "C", value=c)
    lseries = Part("pyspice", "L", value=l)

    # ========================== connections =========================
    rseries[1] += cparallel[1]
    rseries[2] += cparallel[2]
    rseries[1] += lseries[2]
    return Subcircuit(pins=[lseries[1], rseries[2]])
