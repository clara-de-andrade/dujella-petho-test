from mpmath import nint
from math import floor, ceil


def nint_dist(x):
    return abs(x - nint(x))

def prog_bar(i, i_max):
    prog = int(100 * i / i_max)
    return f"[{"#" * floor(prog / 2)}{" " * ceil((100 - prog) / 2)}] {prog}%   "