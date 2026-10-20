#!/usr/bin/python3

from mpmath import mp as ctx
from mpmath import mpf as mpf
from mpmath import nstr, exp, log, power, floor, ceil, findroot

from base import nint_dist, prog_bar
from cache import load_cache, update_cache, purge_cache


def init_cache(prec, m_max):
    print(
        "Inicializando cachê com parâmetros:",
        f"  prec  = {prec}",
        f"  m_max = {m_max}",
        "",
        sep="\n"
    )

    cache = load_cache()
    flag = False

    if prec != cache.get("prec"):
        cache = purge_cache()
        flag = True

    # Precisão dos cálculos com números decimais no mpmath, em bits.
    ctx.prec = prec

    if m_max != cache.get("m_max"):
        cache = purge_cache()
        flag = True

    s_max = cache.get("s_max")
    if not s_max:
        # Cota superior de s. Dada por s <= 4.7e14 * (m ** 5) * (log(m) ** 4).
        s_max = floor(mpf("4.7e14") * (m_max ** 5) * (log(m_max) ** 4))    
        flag = True

    k_max = cache.get("k_max")
    if not k_max:
        # Cota superior de k. Dada por k <= log(s).
        k_max = int(floor(log(s_max)))
        flag = True

    n_max = cache.get("n_max")
    if not n_max:
        # Cota superior de n. Dada por n <= m * s + 2.
        n_max = m_max * s_max + 2
        flag = True

    if True:
        cache = update_cache({
            "prec": prec,
            "m_max": m_max,

            "n_max": n_max,
            "s_max": s_max,
            "k_max": k_max,
        })
    
    return cache


def init_kbonacci(cache):
    m_max = int(cache["m_max"])
    k_max = int(cache["k_max"])
    
    print(
        f"Pré-computando F(k, m + 1) para",
        f"\t3 <= m <= {m_max}",
        f"\t3 <= k <= {k_max}",
        sep="\n")

    kbonacci = cache.get("kbonacci", {})
    if not kbonacci:
        for k in range(3, k_max + 1):
            prev = [0] * (k - 2) + [1] * 2
            for m in range(3, m_max + 2):
                term = sum(prev)
                kbonacci[m, k] = term
                prev = prev[1:] + [term]

            print(prog_bar(k, k_max), end="\r")
        print()

        cache = update_cache({
            "kbonacci": kbonacci
        })

    print(
        "PRONTO!",
        "",
        sep="\n")

    return cache


def init_alpha_g(cache):
    k_max = int(cache["k_max"])
    prec = cache["prec"]

    print(
        f"Pré-computando alpha_k e g_k para",
        f"\t3 <= k <= {k_max}",
        f"com precisão de {prec} bits",
        sep="\n")

    alpha = cache.get("alpha", {})
    g = cache.get("g", {})

    if not alpha or not g:
        for k in range(3, k_max + 1):
            p = lambda x: x**k - sum(x**i for i in range(k))
            alpha[k] = findroot(p, 2)

            g[k] = (alpha[k] - 1) / (2 + (k + 1) * (alpha[k] - 2))

            print(prog_bar(k, k_max), end="\r")
        print()

        cache = update_cache({
            "alpha": alpha,
            "g": g
        })

    print(
        "PRONTO!",
        "",
        sep="\n")

    return cache


def init_gamma_mu(cache):
    m_max = int(cache["m_max"])
    k_max = int(cache["k_max"])
    
    kbonacci = cache["kbonacci"]
    alpha    = cache["alpha"]
    g        = cache["g"]

    prec = cache["prec"]

    print(
        f"Pré-computando gamma(k, m) e mu(k, m) para",
        f"\t3 <= m <= {m_max}",
        f"\t3 <= k <= min(m, {k_max})",
        f"com precisão de {prec} bits",
        sep="\n")

    gamma = cache.get("gamma", {})
    mu = cache.get("mu", {})

    if not gamma or not mu:
        for m in range(3, m_max + 1):
            for k in range(3, min(m, k_max) + 1):
                log_kbonacci = log(kbonacci[m + 1, k])
                log_alpha = log(alpha[k])
                log_g = log(g[k])

                gamma[m, k] = log_alpha / log_kbonacci
                mu[m, k] = (log_alpha - log_g) / log_kbonacci

            print(prog_bar(m, m_max), end="\r")
        print()
        
        cache = update_cache({
            "gamma": gamma,
            "mu": mu
        })

    print(
        "PRONTO!",
        "",
        sep="\n")

    return cache


def test_q_epsilon(cache):
    m_max = int(cache["m_max"])
    n_max = cache["n_max"]
    k_max = int(cache["k_max"])

    n_max = cache["n_max"]
    gamma = cache["gamma"]
    mu    = cache["mu"]

    prec = cache["prec"]

    print(
        f"Testando q(k, m) > 6 * n e epsilon(k, m) > 0 para",
        f"\t3 <= m <= {m_max}",
        f"\t3 <= k <= min(m, {k_max})",
        f"\tn <= {nstr(n_max, 3)}",
        f"com precisão de {prec} bits",
        sep="\n")

    q = cache.get("q", {})
    epsilon = cache.get("epsilon", {})

    if not q or not epsilon:
        for m in range(3, m_max + 1):
            for k in range(3, min(m, k_max) + 1):
                q_old = 1
                q[m, k] = 0
                eta = gamma[m, k]

                while True:
                    a = floor(eta)
                    eta = 1 / (eta - a)
                    q_old, q[m, k] = q[m, k], a * q[m, k] + q_old
                    if not q[m, k] > 6 * n_max:
                        continue

                    epsilon[m, k] = nint_dist(q[m, k] * mu[m, k]) - n_max * nint_dist(q[m, k] * gamma[m, k])
                    if not epsilon[m, k] > 0:
                        continue

                    break

            print(prog_bar(m, m_max), end="\r")
        print()
        
        cache = update_cache({
            "q": q,
            "epsilon": epsilon
        })

    print(
        "PRONTO!",
        "",
        sep="\n")

    return cache


def init_cotas(cache):
    m_max = cache["m_max"]
    n_max = cache["n_max"]
    k_max = cache["k_max"]
    s_max = cache["s_max"]


    q       = cache["q"]
    epsilon = cache["epsilon"]

    prec = cache["prec"]

    print(
        f"Computando novas cotas para n, k, s",
        f"com precisão de {prec} bits",
        sep="\n")

    log_A = log(mpf("2.01"))
    log_B = log(mpf("1.65"))
    
    n_max_new = -1
    for m in range (3, m_max + 1):
        for k in range(3, min(m, k_max) + 1):
            # Candidato a nova cota superior de n. Dado por
            #   n <= n_max * log(A * q(k, m_t / epsilon(k, m)_t) / log B.
            temp = int(floor(m_max * (log_A + log(q[m, k]) - log(epsilon[m, k])) / log_B))
            if n_max_new < 0 or temp < n_max_new:
                n_max_new = temp

        print(prog_bar(m, m_max), end="\r")
    print()
    
    # Nova cota superior de s. Dada por s <= (n - 2) * (m - 1).
    s_max_new = int(floor((n_max_new - 2) / (3 - 1)))
    # Nova cota superior de k. Dada por k <= log s.
    k_max_new = int(floor(log(s_max_new)))

    if [n_max, k_max, s_max] != [n_max_new, k_max_new, s_max_new]:
        cache = update_cache({
            "n_max": n_max_new,
            "s_max": s_max_new,
            "k_max": k_max_new,
        })

    print(
        "PRONTO!",
        "",
        sep="\n")

    return cache


