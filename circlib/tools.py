from typing import List
import math


def noise_figure_equiv(inputs: List[List[float]]) -> float:
    """
    Calculate the equivalent noise figure for multiple devices chained
    together.

    :param inputs: A 2D list, where the outer dimensions corresponds
        to each device in the chain and the inner dimension has two
        elements, where the first element is the noise figure for that
        device and the second element is the gain.  Both of these
        values must be given in dB.

    :returns: The total equivalent noise figure for all devices (in
              dB).
    """
    n = len(inputs)
    f = [10 ** (elem[0] / 10) for elem in inputs]
    g = [10 ** (elem[1] / 10) for elem in inputs]
    fnet = f[0]
    gnet = g[0]
    for i in range(1, n):
        fnet += (f[i] - 1) / gnet
        gnet += g[i]

    return 10 * math.log(fnet, 10)
