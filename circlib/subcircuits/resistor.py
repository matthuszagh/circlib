"""
"""
from circlib.quantity import (
    Power,
    Resistance,
    Inductance,
    Capacitance,
)
from circlib.eseries import Eseries
from circlib.net import Net
from circlib.kicad.footprint import KicadFootprint, ResistorSMD
from circlib.passive import Passive
from skidl import Part


class Resistor(Passive):
    """
    """

    def __init__(
        self,
        nom: float = None,
        eseries: Eseries = Eseries(series=24, tol=1e-2, limits=(1, 1e7)),
        power_rating: Power = Power(typval=1 / 4),
        footprint: KicadFootprint = ResistorSMD["0402"],
        model: str = None,
    ):
        """
        :param nom: Nominal resistance value.  The actual resistance
            value will be set to the closest E-series value for the
            provided ``eseries`` parameter.  If ``eseries`` is set to
            None, the resistance will be set exactly to the value of
            ``nom``.  ``nom`` can also be set to None if you're only
            using the resistor as a way to direct other subcircuits.
        :param series: The E-series standard.  This can be set to None
            to make the resistance exactly equal to ``nom``.
        :param power_rating: Resistor power rating.
        :param footprint: KiCAD footprint.
        :param model: Spice model.  If left as the default None, a
            single resistance value (i.e. no parasitics) will be used.
            This can be set to any two-port network spice netlist for
            more sophisticated parasitic simulation.
        """
        self.eseries = eseries
        if nom is not None:
            if self.eseries is None:
                self.val = Resistance(nom)
            else:
                self._set_val(nom)
        self.power_rating = power_rating
        self.footprint = footprint
        self.model = model

    def construct(
        self, n1: Net = None, n2: Net = None,
    ):
        """
        :param n1: Net tied to node 1 of the resistor.  If left as
            None, this will have to be set later.
        :param n2: Net tied to node 2 of the resistor.  If left as
            None, this will have to be set later.
        """
        self.cad = Part("Device", "R", footprint=self.footprint.to_skidl())

    def resistance(self, simplify: bool = True) -> Resistance:
        """
        """
        return self.val

    def inductance(self, simplify: bool = True) -> Inductance:
        """
        """
        return Inductance(0)

    def capacitance(self, simplify: bool = True) -> Capacitance:
        """
        """
        return Capacitance(0)

    def _set_val(self, nom: float) -> None:
        """
        Set the resistance value based on the E-series and nominal.
        """
        typval = self.eseries.nearest(nom)
        minval = (1 - self.eseries.tol) * typval
        maxval = (1 + self.eseries.tol) * typval
        self.val = Resistance(typval=typval, minval=minval, maxval=maxval)
