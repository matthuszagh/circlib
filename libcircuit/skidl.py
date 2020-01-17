from skidl import Pin
from typing import List


class Subcircuit:
    def __init__(self, pins: List[Pin]):
        self.pins = pins

    def __getitem__(self, key: int):
        return self.pins[key - 1]
        # return getattr(self.pins, self.pins[key - 1])

    def __setitem__(self, key: int, value):
        self.pins[key - 1] = value
