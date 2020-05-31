from typing import List, Union
from abc import ABC
from functools import partial


# TODO all these classes aren't necessary, since they are mostly
# identical. The only difference generally is the type of value we'd
# like to represent.


class Quantity(ABC):
    """
    A physical quantity.  All values should be given in SI units (e.g.
    a value of 1 for current signifies 1A).

    Quantity uses the datasheet convention of specifying maximum,
    minimum and typical values.  This provides a bit more flexibility
    than the approach of using a nominal value and associated error
    tolerance.

    To support the math operations on Quantity that propagate errors,
    all values (``typval``, ``minval``, ``maxval``) must be internally
    defined.  However, it's only required that you define one of these
    (usually ``typval``, although it technically makes no difference).
    In this case, the other values will be set to the value you
    define.  If you define two values, the definition of the third is
    dependent on the two you define.  If you define ``minval`` and
    ``maxval``, ``typval`` will be placed equidistant between these
    two.  If you define ``typval`` and either ``minval`` or
    ``maxval``, the omitted parameter will be set to ``typval``.
    """

    def __init__(
        self,
        typval: Union[float, complex] = None,
        minval: Union[float, complex] = None,
        maxval: Union[float, complex] = None,
        conditions: List[str] = None,
    ):
        """
        :param typval: Typical value.
        :param minval: Minimum value.
        :param maxval: Maximum value.
        :param conditions: List of test conditions for which
            ``minval``, ``maxval`` and ``typval`` are valid.  In its
            current state, this parameter is purely informational.
            That is, it is unused in the code base.
        """
        args = [typval, minval, maxval]
        none_count = args.count(None)
        if none_count == 3:
            raise ValueError(
                "You must provide at least 1 of typval, minval, and maxval."
            )
        elif none_count == 2:
            non_nulls = [val for val in args if val is not None]
            val = non_nulls[0]
            self.typval = val
            self.minval = val
            self.maxval = val
        elif none_count == 1:
            if typval is None:
                self.typval = (minval + maxval) / 2
                self.minval = minval
                self.maxval = maxval
            elif minval is None:
                self.minval = typval
                self.typval = typval
                self.maxval = maxval
            else:
                self.maxval = typval
                self.typval = typval
                self.minval = minval
        else:
            self.typval = typval
            self.minval = minval
            self.maxval = maxval

        self._conditions = conditions

    def __add__(self, other):
        return Op("+", self, other)

    def __sub__(self, other):
        return Op("-", self, other)

    def __mul__(self, other):
        return Op("*", self, other)

    def __truediv__(self, other):
        return Op("/", self, other)

    def __pow__(self, other):
        return Op("**", self, other)


class Const(Quantity):
    """
    """

    def __init__(self, val: float):
        """
        """
        super().__init__(minval=val, maxval=val, typval=val)


class ComplexQuantity(Quantity, ABC):
    """
    """

    def __init__(
        self,
        minval: complex = None,
        maxval: complex = None,
        typval: complex = None,
        conditions: List[str] = None,
    ):
        """
        """
        super().__init__(minval=minval, maxval=maxval, typval=typval)


class DCQuantity(Quantity, ABC):
    """
    """

    def __init__(
        self,
        minval: float = None,
        maxval: float = None,
        typval: float = None,
        conditions: List[str] = None,
    ):
        """
        """
        super().__init__(minval=minval, maxval=maxval, typval=typval)


class ACQuantity(Quantity, ABC):
    """
    """

    def __init__(
        self,
        minval: float = None,
        maxval: float = None,
        typval: float = None,
        conditions: List[str] = None,
        dc: DCQuantity = None,
    ):
        """
        All AC parameters (``minval``, ``maxval``, ``typval``) are
        interpreted as amplitudes.

        :param dc: DC quantity offset.
        """
        super().__init__(minval=minval, maxval=maxval, typval=typval)
        self._dc = dc

    @property
    def dc(self):
        """
        """
        return self._dc


class DCCurrent(DCQuantity):
    """
    """

    def __init__(
        self,
        minval: float = None,
        maxval: float = None,
        typval: float = None,
        conditions: List[str] = None,
    ):
        """
        """
        super().__init__(
            minval=minval, maxval=maxval, typval=typval, conditions=conditions
        )


class ACCurrent(ACQuantity):
    """
    """

    def __init__(
        self,
        minval: float = None,
        maxval: float = None,
        typval: float = None,
        conditions: List[str] = None,
        dc: DCCurrent = None,
    ):
        """
        """
        super().__init__(
            minval=minval,
            maxval=maxval,
            typval=typval,
            conditions=conditions,
            dc=dc,
        )


class DCVoltage(DCQuantity):
    """
    """

    def __init__(
        self,
        minval: float = None,
        maxval: float = None,
        typval: float = None,
        conditions: List[str] = None,
    ):
        """
        """
        super().__init__(
            minval=minval, maxval=maxval, typval=typval, conditions=conditions
        )


class ACVoltage(ACQuantity):
    """
    """

    def __init__(
        self,
        minval: float = None,
        maxval: float = None,
        typval: float = None,
        conditions: List[str] = None,
        dc: DCVoltage = None,
    ):
        """
        """
        super().__init__(
            minval=minval,
            maxval=maxval,
            typval=typval,
            conditions=conditions,
            dc=dc,
        )


class Power(DCQuantity):
    """
    """

    def __init__(
        self,
        minval: float = None,
        maxval: float = None,
        typval: float = None,
        conditions: List[str] = None,
    ):
        """
        """
        super().__init__(
            minval=minval, maxval=maxval, typval=typval, conditions=conditions
        )


class Resistance(Quantity):
    """
    """

    def __init__(
        self,
        minval: float = None,
        maxval: float = None,
        typval: float = None,
        conditions: List[str] = None,
    ):
        """
        """
        super().__init__(
            minval=minval, maxval=maxval, typval=typval, conditions=conditions
        )


class Capacitance(Quantity):
    """
    """

    def __init__(
        self,
        minval: float = None,
        maxval: float = None,
        typval: float = None,
        conditions: List[str] = None,
    ):
        """
        """
        super().__init__(
            minval=minval, maxval=maxval, typval=typval, conditions=conditions
        )


class Inductance(Quantity):
    """
    """

    def __init__(
        self,
        minval: float = None,
        maxval: float = None,
        typval: float = None,
        conditions: List[str] = None,
    ):
        """
        """
        super().__init__(
            minval=minval, maxval=maxval, typval=typval, conditions=conditions
        )


# class Impedance(ComplexQuantity):
#     """
#     """

#     def __init__(
#         self,
#         minval: complex = None,
#         maxval: complex = None,
#         typval: complex = None,
#         conditions: List[str] = None,
#     ):
#         """
#         """
#         super().__init__(
#             minval=minval, maxval=maxval, typval=typval, conditions=conditions
#         )


class Temperature(DCQuantity):
    """
    """

    def __init__(
        self,
        minval: float = None,
        maxval: float = None,
        typval: float = None,
        conditions: List[str] = None,
    ):
        """
        All values are assumed to be in Kelvin.  Use a conversion
        function if you need Farenheit or Celcius.
        """
        super().__init__(
            minval=minval, maxval=maxval, typval=typval, conditions=conditions
        )


class Op:
    """
    """

    def __init__(self, op, a, b):
        """
        """
        self.op = op
        self.a = a
        self.b = b

    def __add__(self, other):
        return Op("+", self, other)

    def __sub__(self, other):
        return Op("-", self, other)

    def __mul__(self, other):
        return Op("*", self, other)

    def __truediv__(self, other):
        return Op("/", self, other)

    def __pow__(self, other):
        return Op("**", self, other)


def _vars(expr: Op, indep_vars: List = []) -> List[Quantity]:
    """
    Return all variables in ``expr``.
    """
    a = expr.a
    b = expr.b

    for val in a, b:
        if issubclass(type(val), Quantity):
            if val not in indep_vars:
                indep_vars.append(val)
        else:
            indep_vars = _vars(val, indep_vars)

    return indep_vars


def _func_expr(expr: Op, indep_vars: List, var_names: List) -> str:
    """
    """
    op = expr.op
    a = expr.a
    b = expr.b

    pre_str = "("
    post_str = ")"
    if issubclass(type(a), Quantity):
        i = indep_vars.index(a)
        pre_str += var_names[i]
    else:
        pre_str += _func_expr(a, indep_vars, var_names)

    pre_str += op

    if issubclass(type(b), Quantity):
        i = indep_vars.index(b)
        pre_str += var_names[i]
    else:
        pre_str += _func_expr(b, indep_vars, var_names)

    return pre_str + post_str


def _func(expr: Op, indep_vars: List):
    """
    """
    func_name = "_expr_func"
    arg_names = ["x{}".format(i) for i, _ in enumerate(indep_vars)]
    args = ", ".join(arg_names)
    func_expr = _func_expr(expr, indep_vars, arg_names)
    exec(
        """
def {}({}):
    return {}
    """.format(
            func_name, args, func_expr
        )
    )
    return locals()[func_name]


def _func_partials(func, indep_vars: List):
    """
    """
    partials = []
    for j in range(len(indep_vars)):
        dct = {}
        for i, var in enumerate(indep_vars):
            if not i == j:
                dct["x{}".format(i)] = var.typval
        partials.append(partial(func, **dct))

    return partials


def evaluate(expr: Op):
    """
    """
    indep_vars = _vars(expr)
    func = _func(expr, indep_vars)
    partials = _func_partials(func, indep_vars)
    derivs = []
    i = 0
    for part, var in zip(partials, indep_vars):
        x0 = var.typval
        dx = var.typval / 100
        dy = part(**{"x{}".format(i): x0 + dx}) - part(**{"x{}".format(i): x0})
        derivs.append(dy)
        i += 1

    min_dict = {}
    max_dict = {}
    typ_dict = {}
    i = 0
    for deriv, var in zip(derivs, indep_vars):
        if deriv >= 0:
            min_dict["x{}".format(i)] = var.minval
            max_dict["x{}".format(i)] = var.maxval
        else:
            min_dict["x{}".format(i)] = var.maxval
            max_dict["x{}".format(i)] = var.minval
        typ_dict["x{}".format(i)] = var.typval
        i += 1

    minval = func(**min_dict)
    maxval = func(**max_dict)
    typval = func(**typ_dict)

    return Quantity(minval=minval, maxval=maxval, typval=typval)
