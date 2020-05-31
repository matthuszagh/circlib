# from skidl import (
#     Part,
#     Interface,
#     Pin,
#     lib_search_paths,
#     footprint_search_paths,
# )


# class L7980(Interface):
#     """
#     """

#     def __init__(self):
#         """
#         """


# class L7980:
#     """
#     TODO: iout_min/max should be calculated based on downstream
#     devices.  We should be able to use the output net in some way to
#     do this.
#     """

#     def __init__(
#         self,
#         vinp: Pin,
#         voutp: Pin,
#         enp: Pin,
#         syncp: Pin,
#         gndp: Pin,
#         vout_target: float,
#         fsw_target: float,
#         iout_min: float,
#         iout_max: float,
#         vin_nominal: float,
#         vin_pp: float,
#         diode: Diode = None,
#     ):
#         """
#         :param vin_pp: Peak-to-peak input voltage.
#         """
#         # pins
#         self._vinp = vinp
#         self._voutp = voutp
#         self._enp = enp
#         self._syncp = syncp
#         self._gndp = gndp

#         # parts
#         self._l7980 = None
#         self._c_bypass = None
#         self._c_bulk = None
#         self._r_sw = None
#         self._r1 = None
#         self._r2 = None

#         self._vout_target = vout_target
#         self._fsw_target = fsw_target
#         self._fsw = None
#         self._check_fsw()
#         self._vin_nominal = vin_nominal
#         self._vin_pp = vin_pp
#         self._iout_min = iout_min
#         self._iout_max = iout_max
#         self._check_output_current()

#         self._construct()

#     def _check_output_current(self):
#         """
#         """
#         if self._iout_max > 2:
#             raise ValueError(
#                 "L7980 only supports a maximum output current of 2A."
#             )

#     @subcircuit
#     def _construct(self):
#         """
#         """
#         self._l7980 = Part(
#             "fmcw",
#             "l7980",
#             footprint="fmcw_v5:HSOP-8-1EP_3.9x4.9mm_P1.27mm_EP2.3x2.3mm_ThermalVias",
#             mfn="L7980A",
#         )
#         self._connect_ground()
#         self._construct_bypass_cap()
#         self._construct_bulk_cap()
#         self._construct_r_sw()

#     def _connect_ground(self):
#         """
#         """
#         self._l7980[7] += self._l7980[9], self._gndp

#     def _construct_bypass_cap(self):
#         """
#         """
#         self._c_bypass = Part(
#             "Device",
#             "C",
#             footprint="Capacitor_SMD:C_0402_1005Metric",
#             value=100e-9,
#         )
#         self._l7980[8] += self._c_bypass[1]
#         self._c_bypass[2] += self._gndp

#     def _construct_bulk_cap(self):
#         """
#         """
#         self._c_bulk = Part(
#             "Device",
#             "C",
#             footprint="Capacitor_SMD:C_0805_2012Metric",
#             value=10e-6,
#         )
#         self._l7980[8] += self._c_bulk[1]
#         self._c_bulk[2] += self._gndp

#     def _construct_r_sw(self):
#         """
#         """
#         if self._fsw_target < 250e3 or self._fsw_target > 1e6:
#             raise ValueError(
#                 "Switching frequency must be in the range from 250kHz to 1MHz."
#             )

#         rsw_val = 28.5e9 / (self._fsw_target - 250e3) - 3.23e3
#         std_val = eseries_val(rsw_val)
#         self._fsw = 28.5e9 / (std_val + 3.23e3) + 250e3
#         print(
#             "Actual switching frequency set to {:.2f}kHz using standard resistor value.".format(
#                 self._fsw / 1e3
#             )
#         )
#         self._r_sw = Part(
#             "Device",
#             "R",
#             footprint="Resistor_SMD:R_0402_1005Metric",
#             value=std_val,
#         )
#         self._l7980[6] += self._r_sw[1]
#         self._r_sw[2] += self._gndp

#     def _min_inductance(self):
#         """
#         """
