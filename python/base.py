from mpmath import nint
from math import floor, ceil


def nint_dist(x):
    return abs(x - nint(x))

def qconv(eta, t):
    u, v = 1, 0

    for _ in range(t):
        a = floor(eta)
        u, v = v, a * v + u
        eta = 1 / (eta - a)

    return v

def prog_bar(i, i_max):
    prog = int(100 * i / i_max)
    return f"[{"#" * floor(prog / 2)}{" " * ceil((100 - prog) / 2)}] {prog}%   "