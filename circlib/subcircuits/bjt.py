"""
"""

from circlib.quantity import DCVoltage, DCCurrent, Temperature
from circlib.const import room_temperature
from circlib.unit import c_to_k
from circlib.net import Net


class NPN:
    """
    """

    def __init__(self, model: str = None, spice: bool = True):
        """
        :param model: Name of the spice model to use for simulating
            the transistor's behavior.  If None is used, a default
            model without parasitics is used.
        :param spice: Whether to use spice to simulate the NPN
            transistor behavior.  If set to False, approximate
            equations will be used instead.  This will be faster but
            is less accurate.
        """

    def vbe(
        self,
        ic: DCCurrent,
        temp: Temperature = Temperature(
            c_to_k(20), c_to_k(80), room_temperature()
        ),
    ) -> DCVoltage:
        """
        Base-emitter voltage drop.
        """

    def construct(self, nc: Net, nb: Net, ne: Net, ns: Net = None):
        """
        """
