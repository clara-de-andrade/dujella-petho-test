from mpmath import nint


def nint_dist(x):
    return abs(x - nint(x))

def prog_bar(i, i_max):
    prog = int(100 * i / i_max)
    return f"|{"=" * (prog // 5)}{"." * ((100 - prog) // 5)}| {prog}%   "