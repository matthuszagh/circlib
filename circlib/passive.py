"""
"""

from __future__ import annotations
from typing import Union, List
from abc import ABC, abstractmethod
from circlib.quantity import (
    Resistance,
    Inductance,
    Capacitance,
    evaluate,
    Const,
    Op,
)


class Passive(ABC):
    """
    Abstract base class for a passive component.  All basic passive
    components (resistors, capacitors, and inductors) inherit from
    this class.
    """

    @abstractmethod
    def resistance(self, simplify: bool = True) -> Union[Resistance, Op]:
        """
        """

    @abstractmethod
    def inductance(self, simplify: bool = True) -> Union[Inductance, Op]:
        """
        """

    @abstractmethod
    def capacitance(self, simplify: bool = True) -> Union[Capacitance, Op]:
        """
        """


class PassiveNetwork(ABC):
    """
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
        self.construct()

    @classmethod
    def from_passives(cls, passives: List[Passive]):
        """
        Build a parallel or series network from a list of passives.
        The order of passives in the network will match that of the
        order in the provided list.
        """
        if len(passives) < 2:
            raise ValueError(
                "List of passives must contain at least 2 elements."
            )

        first = passives.pop(0)
        second = passives.pop(0)
        network = cls(e1=first, e2=second)

        while len(passives) != 0:
            nxt = passives.pop(0)
            network = cls(e1=network, e2=nxt)

        return network

    @abstractmethod
    def construct(self):
        """
        """

    def resistor_count(self) -> int:
        """
        Number of resistors in the network.
        """
        count = 0
        for subent in [self.e1, self.e2]:
            if issubclass(subent, Passive) and subent.resistance().typval != 0:
                count += 1
            else:
                count += subent.resistor_count()

        return count

    def inductor_count(self) -> int:
        """
        Number of inductors in the network.
        """
        count = 0
        for subent in [self.e1, self.e2]:
            if issubclass(subent, Passive) and subent.inductance().typval != 0:
                count += 1
            else:
                count += subent.inductance_count()

        return count

    def capacitor_count(self) -> int:
        """
        Number of capacitors in the network.
        """
        count = 0
        for subent in [self.e1, self.e2]:
            if (
                issubclass(subent, Passive)
                and subent.capacitance().typval != 0
            ):
                count += 1
            else:
                count += subent.capacitance_count()

        return count

    @abstractmethod
    def resistance(self, simplify: bool = True) -> Union[Resistance, Op]:
        """
        """

    @abstractmethod
    def inductance(self, simplify: bool = True) -> Union[Inductance, Op]:
        """
        """

    @abstractmethod
    def capacitance(self, simplify: bool = True) -> Union[Capacitance, Op]:
        """
        """


class ParallelNetwork(PassiveNetwork):
    """
    A network of two entities in parallel.  An entity can either be a
    Passive or PassiveNetwork.  All Passives must have their
    ``construct`` member functions called prior to being added to this
    network.
    """

    def construct(self):
        """
        """
        for i in range(1, 3):
            self.e1.cad[i] += self.e2.cad[i]

    def resistance(self, simplify: bool = True) -> Union[Resistance, Op]:
        """
        """
        if (
            self.e1.resistance().typval == 0
            or self.e2.resistance().typval == 0
        ):
            return Resistance(0)

        expr = (
            Const(1) / self.e1.resistance(simplify=False)
            + Const(1) / self.e2.resistance(simplify=False)
        ) ** Const(-1)
        if simplify:
            return evaluate(expr)
        return expr

    def inductance(self, simplify: bool = True) -> Union[Inductance, Op]:
        """
        """
        if (
            self.e1.inductance().typval == 0
            or self.e2.inductance().typval == 0
        ):
            return Inductance(0)

        expr = (
            Const(1) / self.e1.inductance(simplify=False)
            + Const(1) / self.e2.inductance(simplify=False)
        ) ** Const(-1)
        if simplify:
            return evaluate(expr)
        return expr

    def capacitance(self, simplify: bool = True) -> Union[Capacitance, Op]:
        """
        """
        expr = self.e1.capacitance(simplify=False) + self.e2.capacitance(
            simplify=False
        )
        if simplify:
            return evaluate(expr)
        return expr


class SeriesNetwork(PassiveNetwork):
    """
    A network of two entities in series.  An entity can either be a
    Passive or PassiveNetwork.
    """
