from typing import Tuple, List, Callable, Optional
from bisect import bisect_left
import numpy as np


def normalize(val: float, multiplier: float = 1) -> Tuple[float, float]:
    """
    Find two numbers, a and b, such that a*b=val, a is in the range
    [1,10), and b%10==0.

    :param val: Value to normalize to the range [1,10).  Must be > 0.
    :param multiplier: Value multiplier.  The value encoded in this
        function is always equal to ``val``*``multiplier``.

    :returns: (a, b) in the description, or the normalized value and
              its multiplier.
    """
    if val >= 1 and val < 10:
        return (val, multiplier)

    if val >= 10:
        val /= 10
        multiplier *= 10
        return normalize(val, multiplier)

    if val < 1:
        val *= 10
        multiplier /= 10
        return normalize(val, multiplier)


def _nearest_in_list(val: float, lst: List[float]) -> float:
    """
    """
    ins_pos = bisect_left(lst, val)
    if ins_pos == 0:
        return lst[0]
    if ins_pos == len(lst):
        return lst[-1]
    lower = lst[ins_pos - 1]
    higher = lst[ins_pos]
    if higher - val < val - lower:
        return higher
    else:
        return lower


def _nearest_less_in_list(val: float, lst: List[float]) -> float:
    """
    """
    ins_pos = bisect_left(lst, val)
    if ins_pos == 0 and val != lst[0]:
        raise ValueError("Supplied value is less than all list values.")
    if ins_pos == len(lst):
        return lst[-1]
    if val == lst[ins_pos]:
        return val
    return lst[ins_pos - 1]


def _nearest_greater_in_list(val: float, lst: List[float]) -> float:
    """
    """
    ins_pos = bisect_left(lst, val)
    if ins_pos == len(lst):
        raise ValueError("Supplied value is greater than all list values.")
    return lst[ins_pos]


class Eseries:
    """
    Standard values for use with resistors, capacitors, inductors and
    zener diodes.
    """

    def __init__(
        self,
        series: int = 24,
        tol: float = None,
        limits: Tuple[float, float] = (1, 10e6),
    ):
        """
        :param series: E-series to use.  This can be one of [3, 6, 12,
            24, 48, 96, or 192].
        :param tol: Value tolerance.  This can be one of [40, 20, 10,
            5, 2, 1, 0.5, 0.25, 0.1]e-2.  If left as the default value
            of None, the tolerance will be set to the default
            tolerance for the series.  Note that this will result in
            excessively large tolerances for most cases.
        :param limits: Lower and upper value bounds, inclusive.
        """
        self._set_series(series)
        self._set_tol(tol)
        self._set_limits(limits)

    def nearest(self, val: float) -> float:
        """
        Find the nearest E-series standard value for a given value.

        :param val: Target value.
        """
        if val < self.limits[0]:
            return self.limits[0]
        if val > self.limits[1]:
            return self.limits[1]

        return self._nearest_impl(val, _nearest_in_list)

    def nearest_less(self, val: float) -> float:
        """
        Find the nearest E-series standard value strictly less than
        or equal to the supplied value.
        """
        if val < self.limits[0]:
            raise ValueError(
                "Value requested is less than smallest acceptable value."
            )
        if val > self.limits[1]:
            return self.limits[1]

        return self._nearest_impl(val, _nearest_less_in_list)

    def nearest_greater(self, val: float) -> float:
        """
        Find the nearest E-series standard value strictly greater than
        or equal to the supplied value.
        """
        if val < self.limits[0]:
            return self.limits[0]
        if val > self.limits[1]:
            raise ValueError(
                "Value requested is greater than largest acceptable value."
            )

        return self._nearest_impl(val, _nearest_greater_in_list)

    def next_lower(self, val: float) -> float:
        """
        Find the next E-series value below the provided one.

        :param val: E-series value for which you'd like the next lower
            value.
        """
        if val < self.limits[0] or val > self.limits[1]:
            raise ValueError("Provided value outside limits.")

        norm = normalize(val)
        coeff = round(norm[0], 2)
        mag = round(norm[1], 2)
        if coeff not in self.stdvals:
            raise ValueError("Value is not a standard E-series value.")

        if self.stdvals.index(coeff) == 0:
            mag /= 10
            coeff = round(self.stdvals[-1], 2)
            return round(coeff * mag, 2)

        idx = self.stdvals.index(coeff)
        return round(mag * self.stdvals[idx - 1], 2)

    def next_higher(self, val: float) -> float:
        """
        Find the next E-series value above the provided one.

        :param val: E-series value for which you'd like the next higher
            value.
        """
        if val < self.limits[0] or val > self.limits[1]:
            raise ValueError("Provided value outside limits.")

        norm = normalize(val)
        coeff = round(norm[0], 2)
        mag = round(norm[1], 2)
        if coeff not in self.stdvals:
            raise ValueError("Value is not a standard E-series value.")

        if self.stdvals.index(coeff) == len(self.stdvals) - 1:
            mag *= 10
            coeff = round(self.stdvals[0], 2)
            return round(coeff * mag, 2)

        idx = self.stdvals.index(coeff)
        return round(mag * self.stdvals[idx + 1], 2)

    def _nearest_impl(
        self, val: float, coeff_func: Callable[[float, List[float]], float]
    ) -> float:
        """
        """
        coeff, mul = normalize(val)
        stdcoeff = coeff_func(coeff, self.stdvals)
        return round(stdcoeff * mul, 2)

    def _set_series(self, series: int) -> None:
        """
        """
        valid_series = self._valid_series()
        if series not in valid_series:
            raise ValueError(
                "Invalid series value specified: {}".format(series)
            )
        self.series = series

        if series <= 24:
            stdvals = [round(10 ** (n / series), 1) for n in range(series)]
            # account for official value deviations from calculated value.
            stdvals = [2.7 if x == 2.6 else x for x in stdvals]
            stdvals = [3.0 if x == 2.9 else x for x in stdvals]
            stdvals = [3.3 if x == 3.2 else x for x in stdvals]
            stdvals = [3.6 if x == 3.5 else x for x in stdvals]
            stdvals = [3.9 if x == 3.8 else x for x in stdvals]
            stdvals = [4.3 if x == 4.2 else x for x in stdvals]
            stdvals = [4.7 if x == 4.6 else x for x in stdvals]
            stdvals = [8.2 if x == 8.3 else x for x in stdvals]
        else:
            stdvals = [round(10 ** (n / series), 2) for n in range(series)]

        if series == 192:
            stdvals = [9.20 if x == 9.19 else x for x in stdvals]

        self.stdvals = stdvals

    def _valid_series(self) -> List[int]:
        """
        """
        return [3, 6, 12, 24, 48, 96, 192]

    def _set_tol(self, tol: Optional[float]) -> None:
        """
        """
        valid_tol = self._valid_tol()
        valid_series = self._valid_series()
        if tol is None:
            self.tol = valid_tol[valid_series.index(self.series)]
        else:
            if tol not in valid_tol:
                raise ValueError(
                    "Invalid tolerance value specified {}".format(tol)
                )
            if valid_tol.index(tol) < valid_series.index(self.series):
                raise ValueError(
                    "Tolerance for series {} must be at most {:.2f}%".format(
                        self.series,
                        100 * valid_tol[valid_series.index(self.series)],
                    )
                )
            self.tol = tol

    def _valid_tol(self) -> List[float]:
        """
        """
        return np.multiply(
            [40, 20, 10, 5, 2, 1, 0.5, 0.25, 0.1], 1e-2
        ).tolist()

    def _set_limits(self, limits: Tuple[float, float]) -> None:
        """
        """
        self.limits = (self._lower_lim(limits[0]), self._upper_lim(limits[1]))

    def _lower_lim(self, val: float) -> float:
        """
        Lower value limit.  Unlike the user-provided argument, this is
        guaranteed to be a standard value.

        :param val: Non-standard lower limit.
        """
        coeff, mul = normalize(val)
        stdcoeff = _nearest_greater_in_list(coeff, self.stdvals)
        return round(stdcoeff * mul, 2)

    def _upper_lim(self, val: float) -> float:
        """
        Upper value limit.  Unlike the user-provided argument, this is
        guaranteed to be a standard value.

        :param val: Non-standard lower limit.
        """
        coeff, mul = normalize(val)
        stdcoeff = _nearest_less_in_list(coeff, self.stdvals)
        return round(stdcoeff * mul, 2)
