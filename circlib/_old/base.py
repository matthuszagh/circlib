from abc import ABC, abstractmethod
from libcircuit.skidl import Subcircuit


class BaseCircuit(ABC):
    @abstractmethod
    def cad(self) -> Subcircuit:
        pass

    @abstractmethod
    def spice(self, use_parasitics: bool = True) -> Subcircuit:
        pass

    @abstractmethod
    def _connect_components(self):
        pass
