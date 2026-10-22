#!/usr/bin/python3

import pickle

from mpmath import mp as ctx
from mpmath import mpf as mpf
from mpmath import log, floor, sign, nstr

from base import nint_dist, prog_bar

###############################################################################

m_max = 0
s_max = 0
n_max = 0
k_max = 0

m_min = 3
k_min = 3

def init_cache(prec_init, m_max_init):
    print(
        "Inicializando cachê com parâmetros:",
        f"\tprec  = {prec_init}",
        f"\tm_max = {m_max_init}",
        "",
        sep="\n"
    )

    global m_max
    m_max = m_max_init

    # Precisão dos cálculos com números decimais no mpmath, em bits.
    ctx.prec = prec_init

    # Cota superior de s. Dada por s <= 4.7e14 * (m ** 5) * (log(m) ** 4).
    global s_max
    s_max = floor(mpf("4.7e14") * (m_max ** 5) * (log(m_max) ** 4))    

    # Cota superior de n. Dada por n <= m * s + 2.
    global n_max
    n_max = m_max * s_max + 2

    # Cota superior de k. Dada por k <= log(s).
    global k_max
    k_max = int(floor(log(s_max)))

    print(
        "Cotas iniciais:",
        f"\tn <= {nstr(n_max, 2)}",
        f"\ts <= {nstr(s_max, 2)}",
        f"\tk <= {k_max}",
        "",
        sep="\n"
    )

###############################################################################

F = {}

def init_F():
    global m_max, k_max
    
    print(
        "Pré-computando F(k, m + 1) para",
        f"\t3 <= m <= {m_max}",
        f"\t3 <= k <= {k_max}",
        sep="\n"
    )

    global F

    for k in range(k_min, k_max + 1):
        prev = [0] * (k - 1) + [1] * 1

        F[0, k] = 0
        F[1, k] = 1

        for m in range(2, m_max + 2):
            F[m, k] = sum(prev)
            prev = prev[1:] + [F[m, k]]

        print(prog_bar(k, k_max), end="\r")

    print(
        "",
        "PRONTO!",
        "",
        sep="\n"
    )

###############################################################################

alpha = {}
g = {}

def init_alpha_g():
    global k_max

    print(
        f"Pré-computando alpha_k e g_k para",
        f"\t3 <= k <= {k_max}",
        f"com precisão de {ctx.prec} bits",
        sep="\n")

    global alpha, g

    for k in range(k_min, k_max + 1):
        alpha[k] = mpf(2)
        while True:
            p = (alpha[k] ** k) * (alpha[k] - 2) + 1
            dp = (alpha[k] ** (k - 1)) * ((k + 1) * alpha[k] - 2 * k)
            step = p / dp

            alpha[k] = alpha[k] - step
            if abs(step) < (mpf(2) ** -ctx.prec):
                break

        g[k] = (alpha[k] - 1) / (2 + (k + 1) * (alpha[k] - 2))

        print(prog_bar(k, k_max), end="\r")

    print(
        "",
        "PRONTO!",
        "",
        sep="\n"
    )

###############################################################################

gamma = {}
mu = {}

def init_gamma_mu():
    global k_max, m_max
    
    print(
        "Pré-computando gamma(k, m) e mu(k, m) para",
        f"\t3 <= m <= {m_max}",
        f"\t3 <= k <= min(m, {k_max})",
        f"com precisão de {ctx.prec} bits",
        sep="\n")

    global F, alpha, g, gamma, mu

    for m in range(m_min, m_max + 1):
        for k in range(k_min, min(m, k_max) + 1):
            log_F = log(F[m + 1, k])
            log_alpha = log(alpha[k])
            log_g = log(g[k])

            gamma[m, k] = log_alpha / log_F
            mu[m, k] = (log_alpha - log_g) / log_F

        print(prog_bar(m, m_max), end="\r")

    print(
        "",
        "PRONTO!",
        "",
        sep="\n"
    )

###############################################################################

q = {}
epsilon = {}

def test_q_epsilon():
    global k_max, m_max, n_max

    print(
        "Testando q(k, m) > 6 * n e epsilon(k, m) > 0 para",
        f"\t3 <= m <= {m_max}",
        f"\t3 <= k <= min(m, {k_max})",
        f"\tn <= {nstr(n_max, 2)}",
        f"com precisão de {ctx.prec} bits",
        sep="\n"
    )

    global gamma, mu, q, epsilon

    for m in range(m_min, m_max + 1):
        for k in range(k_min, min(m, k_max) + 1):
            eta = gamma[m, k]
            q_old = 1
            q[m, k] = 0

            while True:
                a = floor(eta)
                eta = 1 / (eta - a)
                q_old, q[m, k] = q[m, k], a * q[m, k] + q_old
                if not q[m, k] > 6 * n_max:
                    continue

                t1 = nint_dist(q[m, k] * mu[m, k])
                t2 = nint_dist(q[m, k] * gamma[m, k])
                epsilon[m, k] = t1 - n_max * t2
                if not epsilon[m, k] > 0:
                    continue

                break

        print(prog_bar(m, m_max), end="\r")
        
    print(
        "",
        "PRONTO!",
        "",
        sep="\n"
    )

###############################################################################

def init_cotas():
    global m_max, k_max

    print(
        "Reduzindo as cotas superiores de n, k, s",
        f"com precisão de {ctx.prec} bits",
        sep="\n"
    )

    log_A = log(mpf("2.01"))
    log_B = log(mpf("1.65"))
    
    # Nova cota superior de n. Dada por
    #   n <= n_max * log(A * q(k, m)_t / epsilon(k, m)_t) / log B.
    global n_max
    for m in range(m_min, m_max + 1):
        for k in range(k_min, min(m, k_max) + 1):
            temp = floor(m_max * (log_A + log(q[m, k]) - log(epsilon[m, k])) / log_B)
            if temp < n_max:
                n_max = int(temp)

        prog_bar(m, m_max)
    
    # Nova cota superior de s. Dada por s <= (n - 2) * (m - 1).
    global s_max
    s_max = int(floor((n_max - 2) / (m_min - 1)))

    # Nova cota superior de k. Dada por k <= log s.
    k_max = int(floor(log(s_max)))

    print(
        "",
        "PRONTO!",
        "",
        sep="\n"
    )

    print(
        "Cotas finais:",
        f"\tn <= {n_max}",
        f"\ts <= {s_max}",
        f"\tk <= {k_max}",
        "",
        sep="\n"
    )


def write_cache():
    cache = {
        "m_max"  : m_max,
        "prec"   : ctx.prec,

        "F"      : F,
        "alpha"  : alpha,
        "g"      : g,
        "gamma"  : gamma,
        "mu"     : mu,

        "q"      : q,
        "epsilon": epsilon,

        "n_max"  : n_max,
        "s_max"  : s_max,
        "k_max"  : m_max,
    }

    print("Escrevendo cachê em 'cache.pkl'")

    with open("cache.pkl", "wb") as f:
        pickle.dump(cache, f)

    print("PRONTO!")
