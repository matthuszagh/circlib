"""
"""

from typing import Union
from abc import ABC, abstractmethod
from circlib.quantity import (
    Resistance,
    Inductance,
    Capacitance,
    Impedance,
    simplify,
    Const,
)


class Passive(ABC):
    """
    Abstract base class for a passive component.  All basic passive
    components (resistors, capacitors, and inductors) inherit from
    this class.
    """

    @abstractmethod
    def resistance(self) -> Resistance:
        """
        """

    @abstractmethod
    def inductance(self) -> Inductance:
        """
        """

    @abstractmethod
    def capacitance(self) -> Capacitance:
        """
        """

    @abstractmethod
    def impedance(self) -> Impedance:
        """
        """


class PassiveNetwork(ABC):
    """
    """


class ParallelNetwork(PassiveNetwork):
    """
    A network of two entities in parallel.  An entity can either be a
    Passive or PassiveNetwork.  All Passives must have their
    ``construct`` member functions called prior to being added to this
    network.
    """

    def __init__(
        self,
        e1: Union[Passive, PassiveNetwork],
        e2: Union[Passive, PassiveNetwork],
    ):
        """
        """
        self.e1 = e1
        self.e2 = e2

    def construct(self):
        """
        """
        for i in range(1, 3):
            self.e1.cad[i] += self.e2.cad[i]

    def resistance(self) -> Resistance:
        """
        """
        if (
            self.e1.resistance().typval == 0
            or self.e2.resistance().typval == 0
        ):
            return Resistance(0)

        return simplify(
            (Const(1) / self.e1.resistance() + Const(1) / self.e2.resistance())
            ** Const(-1)
        )

    def inductance(self) -> Inductance:
        """
        """
        if (
            self.e1.inductance().typval == 0
            or self.e2.inductance().typval == 0
        ):
            return Inductance(0)

        return simplify(
            (Const(1) / self.e1.inductance() + Const(1) / self.e2.inductance())
            ** Const(-1)
        )

    def capacitance(self) -> Capacitance:
        """
        """
        return simplify(self.e1.capacitance() + self.e2.capacitance())

    def impedance(self) -> Impedance:
        """
        """
        if (
            abs(self.e1.impedance().typval) == 0
            or abs(self.e2.impedance().typval) == 0
        ):
            return Impedance(complex(0, 0))

        return simplify(
            (Const(1) / self.e1.impedance() + Const(1) / self.e2.impedance())
            ** Const(-1)
        )


class SeriesNetwork(PassiveNetwork):
    """
    A network of two entities in series.  An entity can either be a
    Passive or PassiveNetwork.
    """

    def __init__(
        self,
        e1: Union[Passive, PassiveNetwork],
        e2: Union[Passive, PassiveNetwork],
    ):
        """
        """
        self.e1 = e1
        self.e2 = e2
        e1.cad[2] += e2.cad[1]
