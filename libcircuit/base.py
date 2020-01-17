from abc import ABC, abstractmethod


class BaseCircuit(ABC):
    @abstractmethod
    def cad(self):
        pass

    @abstractmethod
    def spice(self, use_parasitics: bool = True):
        pass

    @abstractmethod
    def _connect_components(self):
        pass
