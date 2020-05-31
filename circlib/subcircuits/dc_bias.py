"""
"""

from typing import List, Union
from abc import ABC, abstractmethod
from circlib.quantity import (
    Const,
    Quantity,
    DCVoltage,
    Resistance,
    Power,
    DCCurrent,
    evaluate,
    Op,
)
from circlib.net import Net
from circlib.eseries import Eseries, normalize
from circlib.passive import Passive, ParallelNetwork
from circlib.subcircuits.resistor import Resistor


class DCBias(ABC):
    """
    """

    @abstractmethod
    def power(self) -> Power:
        """
        Power dissipation of the dc bias network.
        """

    @abstractmethod
    def construct(self, vp: Net, vn: Net, vout: Net):
        """
        """

    @abstractmethod
    def set_module(self, mod_name: str, new_module):
        """
        """


def optimal_resistors(num: int, target: float, eseries) -> List[float]:
    """
    Find a set of optimal parallel resistor values.  This uses the
    brute-force method of just trying every value in the series with
    the same magnitude.  Although this implementation could clearly be
    improved, it appears to be fast enough even for the 192 series and
    3 resistors.
    """
    if num < 2 or num > 3:
        raise ValueError("Number of parallel resistors can only be 2 or 3.")

    if target > eseries.limits[1]:
        raise ValueError("Target is too large for chosen E-series limits.")

    mag = normalize(target)[1]
    values = [val * mag for val in eseries.stdvals]
    best_vals = [None for _ in range(num)]
    error_fun = lambda x: abs(x - target)

    if num == 2:
        parallel_res_fun = lambda r1, r2: (1 / r1 + 1 / r2) ** -1
        best_error = target
        for val1 in values:
            for val2 in values:
                res = parallel_res_fun(val1, val2)
                error = error_fun(res)
                if error < best_error:
                    best_error = error
                    best_vals[0] = val1
                    best_vals[1] = val2

    if num == 3:
        parallel_res_fun = lambda r1, r2, r3: (1 / r1 + 1 / r2 + 1 / r3) ** -1
        best_error = target
        for val1 in values:
            for val2 in values:
                for val3 in values:
                    res = parallel_res_fun(val1, val2, val3)
                    error = error_fun(res)
                    if error < best_error:
                        best_error = error
                        best_vals[0] = val1
                        best_vals[1] = val2
                        best_vals[2] = val3

    return best_vals


class ResistiveDivider(DCBias):
    """
    A resistive voltage divider.

    The resistor values are chosen to be the maximum values such that
    the output voltage requirement is satisfied for the specified
    current sourcing capability.

    Note that this will not necessarily output the desired voltage in
    a no-load condition.  If you'd like this to be true, set the
    minimum output current to 0.

    TODO: support thermal noise requirements.

    Modules
    -------

        1. ``r_top``: A two-port resistive network acting as the top
           resistor in a resistive voltage divider.  Although this is
           typically just a single resistor, it can be any series or
           parallel two-port network of resistors if desired.

        2. ``r_bot``: A two-port resistive network acting as the
           bottom resistor in a resistive voltage divider.
    """

    def __init__(
        self,
        vp: DCVoltage,
        vn: DCVoltage,
        vout: DCVoltage,
        iout: DCCurrent,
        resistor: Resistor,
        num: int = 2,
        **kwargs,
    ):
        """
        :param vp: Positive voltage supply.
        :param vn: Negative voltage supply.
        :param vout: Desired output voltage.
        :param iout: Output current.
        :param resistor: Resistor to use for each resistor in the
            divider.  If you require different E-series standards or
            resistor footprints for the individual, you must set the
            modules manually.
        :param num: The number of resistors to use.  This must be
            between (inclusive) 1 and 6.  If you supply one of the
            modules (either ``r_top`` or ``r_bot``) this will be
            interpreted as the number of resistors to use for the
            remaining module, with any value greater than 3 being
            interpreted as 3.  If you supply both modules, setting
            this parameter has no effect.  If you set the value to 1
            and do not provide any modules, it will be interpreted as
            2.
        """
        self.vp = vp
        self.vn = vn
        self.vout_set = vout
        self.iout = iout
        self.resistor = resistor
        self.r_top = kwargs.get("r_top")
        self.r_bot = kwargs.get("r_bot")
        if num < 1 or num > 6:
            raise ValueError("Number of resistors must be between 2 and 6.")
        if num == 1:
            if self.r_top is None and self.r_bot is None:
                num = 2
        self.num = num
        self._set_resistors()

    def construct(self, vp: Net, vn: Net, vout: Net):
        """
        """

    def set_module(self, mod_name: str, new_module):
        """
        """

    def resistance_bottom(
        self, simplify: bool = True
    ) -> Union[Resistance, Op]:
        """
        """
        expr = self.r_bot.resistance(simplify=False)
        if simplify:
            return evaluate(expr)
        return expr

    def resistance_top(self, simplify: bool = True) -> Union[Resistance, Op]:
        """
        """
        expr = self.r_top.resistance(simplify=False)
        if simplify:
            return evaluate(expr)
        return expr

    def _set_resistors(self) -> None:
        """
        Set divider resistors.
        """
        if self.r_top is None and self.r_bot is None:
            self._set_top_bottom()
        elif self.r_top is None:
            self._set_top()
        elif self.r_bot is None:
            self._set_bottom()
        else:
            self._resistors_valid()

    def _set_top_bottom(self) -> None:
        """
        Set both resistor values when no modules are provided.
        """
        num_top = self.num // 2
        num_bot = self.num - num_top

        if num_top == num_bot:
            if num_top == 1:
                self._set_top_bottom_singles()

    def _set_top_bottom_singles(self) -> None:
        """
        TODO this is needlessly inefficient.  One way to improve on it
        would be to look at the voltage drop and adjust based on that
        rather than simply decrementing.  That would also provide a
        reasonable path for updating in the case of more than 2
        resistors.
        """
        ratio = self._ratio().typval
        top = self._rt_ubound()
        top = self.resistor.eseries.nearest_less(top)
        self._set_r_top(top)
        bottom = top / ratio
        bottom = self.resistor.eseries.nearest(bottom)
        self._set_r_bot(bottom)
        while not self._resistors_valid():
            top = self.resistor.eseries.next_lower(top)
            self._set_r_top(top)
            bottom = top / ratio
            bottom = self.resistor.eseries.nearest(bottom)
            self._set_r_bot(bottom)

    def _set_r_top(self, nom: float) -> None:
        """
        """
        self.r_top = Resistor(
            nom=nom,
            eseries=self.resistor.eseries,
            power_rating=self.resistor.power_rating,
            footprint=self.resistor.footprint,
        )

    def _set_r_bot(self, nom: float) -> None:
        """
        """
        self.r_bot = Resistor(
            nom=nom,
            eseries=self.resistor.eseries,
            power_rating=self.resistor.power_rating,
            footprint=self.resistor.footprint,
        )

    def _set_top(self) -> None:
        """
        Set the upper resistance network when the bottom is explicitly
        provided.
        """
        if (
            self.r_bot.resistance().maxval * self._ratio().maxval
            > self._rt_ubound()
        ):
            raise ValueError("The bottom resistance provided is too large.")

        target = self.r_bot.resistance().typval * self._ratio().typval
        num = min([self.num, 3])
        values = optimal_resistors(
            num=num, target=target, eseries=self.resistor.eseries
        )
        resistors = [
            Resistor(
                nom=val,
                eseries=self.resistor.eseries,
                power_rating=self.resistor.power_rating,
                footprint=self.resistor.footprint,
            )
            for val in values
        ]
        if num == 1:
            self.r_top = resistors[0]
        else:
            self.r_top = ParallelNetwork.from_passives(resistors)

        if not self._resistors_valid():
            raise RuntimeError(
                "Resistance values do not satisfy the provided constraints. Please fix the bottom module."
            )

    def _set_bottom(self) -> None:
        """
        Set the lower resistance network when the top is explicitly
        provided.
        """
        if self.r_top.resistance().maxval > self._rt_ubound():
            raise ValueError("The top resistance provided is too large.")

        target = self.r_top.resistance().typval / self._ratio().typval
        num = min([self.num, 3])
        values = optimal_resistors(
            num=num, target=target, eseries=self.resistor.eseries
        )
        resistors = [
            Resistor(
                nom=val,
                eseries=self.resistor.eseries,
                power_rating=self.resistor.power_rating,
                footprint=self.resistor.footprint,
            )
            for val in values
        ]
        if num == 1:
            self.r_bot = resistors[0]
        else:
            self.r_bot = ParallelNetwork.from_passives(resistors)

        if not self._resistors_valid():
            raise RuntimeError(
                "Resistance values do not satisfy the provided constraints. Please fix the top module."
            )

    def _resistors_valid(self) -> bool:
        """
        """
        vout = self.vout()
        if (
            vout.maxval > self.vout_set.maxval
            or vout.minval < self.vout_set.minval
        ):
            return False
        return True

    def vout(self, simplify: bool = True) -> Union[Op, DCVoltage]:
        """
        """
        iopen = (self.vp - self.vn) / (
            self.resistance_bottom(simplify=False)
            + self.resistance_top(simplify=False)
        )
        vdrop = self.resistance_top(simplify=False) * (self.iout + iopen)
        expr = self.vp - self.vn - vdrop
        return evaluate(expr) if simplify else expr

    def power(self, simplify: bool = True) -> Union[Op, Power]:
        """
        """
        expr = (self.vp - self.vn) ** Const(2) / (
            self.resistance_bottom(simplify=False)
            + self.resistance_top(simplify=False)
        )
        return evaluate(expr) if simplify else expr

    def _ratio(self, simplify: bool = True) -> Union[Op, Quantity]:
        """
        Top to bottom resistance value ratio.
        """
        expr = (self.vp - self.vn) / self.vout_set - Const(1)
        return evaluate(expr) if simplify else expr

    def _rt_ubound(self) -> float:
        """
        Upper bound estimate for top resistance value.  This does not
        guarantee that all lower values are valid; we still need to
        check each value.  However, all values greater than this value
        are guaranteed to be invalid.
        """
        rt = evaluate(
            (
                (self.vp - self.vn)
                * (Const(1) - Const(1) / (Const(1) + Const(1) / self._ratio()))
                - self.vout_set
            )
            / self.iout
        )
        return rt.maxval
