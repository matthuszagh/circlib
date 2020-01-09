from bisect import bisect_left


def eseries_val(val, series=24):
    """
    Return the value of the nearest E-series standard value. @val is
    the starting desired value and @series is E-series standard to
    use. The E-series can be one of these values: 3, 6, 12, 24, 48,
    96, or 192. It is generally recommended to stick with the default
    value of 24.
    """
    if series not in [3, 6, 12, 24, 48, 96, 192]:
        raise ValueError("Invalid series value specified: {}" % series)

    if val <= 0 or val >= 10:
        raise ValueError(
            "Must provide a value between (non-inclusive) 0 and 10"
        )

    if val <= 24:
        stdvals = [round(10 ** (n / series), 1) for n in range(series)]
        # account for official value deviations from calculated value.
        stdvals = [2.7 if x == 2.6 else x for x in stdvals]
        stdvals = [3.0 if x == 2.9 else x for x in stdvals]
        stdvals = [3.3 if x == 3.2 else x for x in stdvals]
        stdvals = [3.6 if x == 3.5 else x for x in stdvals]
        stdvals = [3.9 if x == 3.8 else x for x in stdvals]
        stdvals = [4.3 if x == 4.2 else x for x in stdvals]
        stdvals = [4.7 if x == 4.6 else x for x in stdvals]
        stdvals = [8.3 if x == 8.2 else x for x in stdvals]
    else:
        stdvals = [round(10 ** (n / series), 2) for n in range(series)]

    mod10 = val % 10
    multiplier = val / mod10

    ins_pos = bisect_left(stdvals, mod10)
    if ins_pos == 0:
        return multiplier * stdvals[0]
    if ins_pos == len(stdvals):
        return multiplier * stdvals[-1]
    lower = stdvals[ins_pos - 1]
    higher = stdvals[ins_pos]
    if higher - mod10 < mod10 - lower:
        return multiplier * higher
    else:
        return multiplier * lower


# TODO this should be improved to support a better tradeoff between
# cost and size.
def resistor_footprint(val, permit_0201=False):
    """
    Returns a recommended footprint for the resistance
    @val. @permit_0201 can be used to allow 0201 footprints to be
    used, which is disabled by default. The returned footprint format
    is a string compatible with KiCad.
    """
    return "Resistors_SMD:R_0402"


# TODO this should be improved to support a better tradeoff between
# cost and size.
def capacitor_footprint(val, permit_0201=False):
    """
    Returns a recommended footprint for the capacitance
    @val. @permit_0201 can be used to allow 0201 footprints to be
    used, which is disabled by default. The returned footprint format
    is a string compatible with KiCad.
    """
    return "Capacitors_SMD:C_0402"
