"""
A collection of amplifier subcircuits.
"""

from abc import ABC, abstractmethod
from circlib.net import Net
from circlib.quantity import ACCurrent, ACVoltage, DCCurrent, DCVoltage
from circlib.calc.bjt import intrinsic_emitter_resistance
from circlib.subcircuits.bjt import NPN


class SingleEndedAmplifier(ABC):
    """
    """

    @abstractmethod
    def construct(self, vin: Net, vout: Net, vp: Net, vn: Net):
        """
        """

    @abstractmethod
    def gain(self):
        """
        """

    @abstractmethod
    def input_impedance(self):
        """
        """

    @abstractmethod
    def output_impedance(self):
        """
        """


class CommonEmitterAmplifier(SingleEndedAmplifier):
    """
    A common emitter amplifier uses an NPN transistor to amplify an AC
    input voltage.  The output voltage is taken from the NPN's
    collector node.

    Modules
    -------

        1. ``dc_bias``: The DC bias sets the DC operating point for
           the AC signal.  The simplest example of this would be a
           resistive divider, although alternatives exist (e.g. a
           matched transistor).

        2. ``dc_block``: The DC block prevents a DC input bias from
           reaching the amplifier input.  The DC bias point is instead
           set by the DC bias module, which allows the bias point to
           operate independently of any upstream devices.  The DC
           block is most commonly just a capacitor.

        3. ``emitter_2p``: The two-port resistive network at the
           amplifying transistor's emitter node.  The simplest
           configurations would be a short circuit for a grounded
           emitter amplifier, or a single resistor for a degenerated
           emitter amplifier.  However, better alternatives exist
           including bypassing schemes.

        4. ``collector_2p``: The two-port resistive network at the
           transistor's collector node between the positive voltage
           supply and the output voltage node.  This is typically just
           a resistor.
    """

    def __init__(
        self,
        vp: DCVoltage,
        vn: DCVoltage,
        vin: ACVoltage,
        vout: ACVoltage,
        iout: DCCurrent,
        npn: NPN,
    ):
        """
        :param vp: Positive supply voltage.
        :param vn: Negative supply voltage.
        :param vin: Input voltage.
        :param vout: Output voltage.
        :param iout: Current draw by all downstream devices.
        :param npn: NPN transistor to use.
        """
        self._rc = None
        self._re = None

    def gain(self):
        """
        Small-signal voltage gain.
        """
        # return -self._rc / (self._re + intrinsic_emitter_resistance())

    def input_impedance(self):
        """
        """

    def output_impedance(self):
        """
        """

    def construct(self, vp: Net, vn: Net, vin: Net, vout: Net):
        """
        """
