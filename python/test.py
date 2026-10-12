
from mpmath import mp as ctx
from mpmath import mpf as mpf
from mpmath import nstr, log, power, floor, findroot

from base import nint_dist, qconv, prog_bar
from cache import load_cache, update_cache, purge_cache


def init_cache(prec, t_min, m_max):
    if prec != load_cache("prec", -1):
        purge_cache()
    if t_min != load_cache("t_min", -1):
        purge_cache()
    if m_max != load_cache("m_max", -1):
        purge_cache()

    # Precisão dos cálculos com números decimais no mpmath, em bits.
    ctx.prec = prec

    # Cota superior de s. Dada por s <= 4.7e14 * (m ** 5) * (log(m) ** 4).
    s_max = floor(mpf("4.7e14") * (m_max ** 5) * (log(m_max) ** 4))    
    # Cota superior de k. Dada por k <= log(s).
    k_max = floor(log(s_max))
    # Cota superior de n. Dada por n <= m * s + 2.
    n_max = m_max * s_max + 2

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

    print(
        "PRONTO!",
        "",
        sep="\n")

    return update_cache({
        "kbonacci": kbonacci
    })


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
            p = lambda x: x**(k + 1) - 2 * x**k + 1
            inf = 2 * (1 - power(2, -k))
            alpha[k] = findroot(p, (inf, 2), solver="anderson")

            g[k] = (alpha[k] - 1) / (2 + (k + 1) * (alpha[k] - 2))

            print(prog_bar(k, k_max), end="\r")
        print()
    print(
        "PRONTO!",
        "",
        sep="\n")

    return update_cache({
        "alpha": alpha,
        "g": g
    })


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

    gamma = load_cache("gamma", {})
    mu = load_cache("mu", {})

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

    print(
        "PRONTO!",
        "",
        sep="\n")

    return update_cache({
        "gamma": gamma,
        "mu": mu
    })


def init_q_epsilon_N(cache):
    m_max = int(cache["m_max"])
    k_max = int(cache["k_max"])
    t_min = int(cache["t_min"])

    n_max = cache["n_max"]
    gamma = cache["gamma"]
    mu    = cache["mu"]

    log_A = log(mpf("2.01"))
    log_B = log(mpf("1.65"))

    print(
        f"Pré-computando q(k, m)_t e epsilon(k, m)_t para",
        f"  3 <= m <= {m_max}",
        f"  3 <= k <= min(m, {k_max})",
        f"  t >= {t_min}",
        sep="\n")

    q = load_cache("q", {})
    epsilon = load_cache("epsilon", {})
    N = load_cache("N", {})

    t = t_min
    if not q or not epsilon or not N:
        successful = False
        while not successful:
            successful = True
            for m in range(3, m_max + 1):
                for k in range(3, min(m, k_max) + 1):
                    q_temp = qconv(gamma[m, k], t)

                    if not q_temp > 6 * n_max:
                        print(f"  q({k}, {m})_{t} == {nstr(q_temp, 10)}")
                        print(f"  q({k}, {m})_{t} <= 6 * n_max")
                        successful = False
                        break

                    q[m, k] = q_temp

                    t1, t2 = nint_dist(q[m, k] * mu[m, k]), nint_dist(q[m, k] * gamma[m, k])
                    if t1 == 0 or t2 == 0:
                        print(f"  Erro de precisão em t = {t}!")
                        successful = False
                        break
        
                    epsilon_temp = t1 - n_max * t2
                    if not epsilon_temp > 0:
                        print(f"  epsilon[{m}, {k}]_{t} == {nstr(epsilon_temp, 10)}")
                        print(f"  epsilon[{m}, {k}]_{t} <= 0")
                        successful = False
                        break

                    epsilon[m, k] = epsilon_temp

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
                # q = {}
                # epsilon = {}
                # N = {}
                # t += 50
        print()
    print(
        "PRONTO!",
        "",
        sep="\n")

    return update_cache({
        "t": t,
        "q": q,
        "epsilon": epsilon,
        "N": N,
    })


def init_cotas(cache):
    N = cache["N"]

    # Nova cota superior de n. Dada por
    #   n <= log(A * q(k, m_t / epsilon(k, m)_t) / log B.
    n_max_out = int(min(N.values()))
    # Nova cota superior de s. Dada por s <= (n - 2) * (m - 1).
    s_max_out = int(floor((n_max_out - 2) / (3 - 1)))
    # Nova cota superior de k. Dada por k <= log s.
    k_max_out = int(floor(log(s_max_out)))

    return update_cache({
        "n_max_out": n_max_out,
        "s_max_out": s_max_out,
        "k_max_out": k_max_out,
    })