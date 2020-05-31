"""
"""


class KicadFootprint:
    """
    """

    def __init__(self, lib: str, fp: str):
        """
        """
        self.lib = lib
        self.fp = fp

    def to_skidl(self) -> str:
        """
        SKiDL-compatible footprint string.
        """
        return self.lib + ":" + self.fp


ResistorSMD = {
    "01005": KicadFootprint("Resistor_SMD", "R_01005_0402Metric"),
    "0201": KicadFootprint("Resistor_SMD", "R_0201_0603Metric"),
    "0402": KicadFootprint("Resistor_SMD", "R_0402_1005Metric"),
    "0603": KicadFootprint("Resistor_SMD", "R_0603_1608Metric"),
    "0805": KicadFootprint("Resistor_SMD", "R_0805_2012Metric"),
    "1210": KicadFootprint("Resistor_SMD", "R_1210_3225Metric"),
}
