
from itertools import product

from mpmath import mp as ctx
from mpmath import mpf as mpf
from mpmath import nstr, exp, log, power, floor, ceil, findroot

from base import nint_dist, qconv, prog_bar
from cache import load_cache, update_cache, purge_cache


def init_cache(prec, t_min, m_max):
    cache = load_cache()
    flag = False

    if prec != cache.get("prec"):
        cache = purge_cache()
        flag = True

    # Precisão dos cálculos com números decimais no mpmath, em bits.
    ctx.prec = prec

    if t_min != cache.get("t_min"):
        cache = purge_cache()
        flag = True
        
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
        return update_cache({
            "prec": prec,
            "t_min": t_min,
            "m_max": m_max,
            "n_max": n_max,
            "s_max": s_max,
            "k_max": k_max,
        })


def init_kbonacci(cache):
    m_max = int(cache["m_max"])
    k_max = int(cache["k_max"])
    
    print(
        f"Pré-computando F(k, m + 1) para",
        f"  3 <= m <= {m_max}",
        f"  3 <= k <= {k_max}",
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

    print(
        f"Pré-computando alpha_k e g_k para",
        f"  3 <= k <= {k_max}",
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

    print(
        f"Pré-computando gamma(k, m) e mu(k, m) para",
        f"  3 <= m <= {m_max}",
        f"  3 <= k <= min(m, {k_max})",
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
    n_max = int(cache["n_max"])
    k_max = int(cache["k_max"])
    t_min = int(cache["t_min"])

    n_max = cache["n_max"]
    gamma = cache["gamma"]
    mu    = cache["mu"]

    log_A = log(mpf("2.01"))
    log_B = log(mpf("1.65"))

    print(
        f"Testando q(k, m)_t > 6 * n e epsilon(k, m)_t > 0 para",
        f"  3 <= m <= {m_max}",
        f"  3 <= k <= min(m, {k_max})",
        f"  t >= {t_min}",
        f"  n <= {n_max}",
        sep="\n")

    q = cache.get("q", {})
    epsilon = cache.get("epsilon", {})
    N = cache.get("N", {})

    t = t_min
    if not q or not epsilon or not N:
        successful = False
        while not successful:
            successful = True
            for m in range(3, m_max + 1):
                for k in range(3, min(m, k_max) + 1):
                    q[m, k] = qconv(gamma[m, k], t)

                    if not q[m, k] > 6 * n_max:
                        print(f"  q({k}, {m})_{t} == {nstr(q[m, k], 10)}")
                        print(f"  q({k}, {m})_{t} <= 6 * n_max")
                        successful = False
                        break

                    t1, t2 = nint_dist(q[m, k] * mu[m, k]), nint_dist(q[m, k] * gamma[m, k])
                    if t1 == 0 or t2 == 0:
                        print(f"  Erro de precisão em t = {t}!")
                        successful = False
                        break
        
                    epsilon[m, k] = t1 - n_max * t2
                    if not epsilon[m, k] > 0:
                        print(f"  epsilon[{m}, {k}]_{t} == {nstr(epsilon[m, k], 10)}")
                        print(f"  epsilon[{m}, {k}]_{t} <= 0")
                        successful = False
                        break

                    log_q = log(q[m, k])
                    log_epsilon = log(epsilon[m, k])
                    N_temp = floor((m_max * (log_A + log_q - log_epsilon) / log_B))

                    if not N_temp > 0:
                        print(f"  Erro de precisão em t = {t}!")
                        successful = False
                        break

                    N[m, k]  = N_temp

                if not successful:
                    break

                print(prog_bar(m, m_max), end="\r")
            
            if successful:
                break
            else:
                exit(1)
        print()
        
        cache = update_cache({
            "t": t,
            "q": q,
            "epsilon": epsilon,
            "N": N
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

    N = cache["N"]
    
    # Nova cota superior de n. Dada por
    #   n <= log(A * q(k, m_t / epsilon(k, m)_t) / log B.
    n_max_new = int(min(N[m, k] for m in range(3, m_max + 1) for k in range(3, min(m, k_max) + 1)))
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

    return cache


