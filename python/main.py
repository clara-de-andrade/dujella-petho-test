#!/usr/bin/python3

from tests import (
    init_cache,
    init_F,
    init_alpha_g,
    init_gamma_mu,
    test_q_epsilon,
    init_cotas
)

if __name__ == "__main__":
    prec  = 4096
    m_max = 1394
    init_cache(prec, m_max)

    init_F()
    init_alpha_g()
    init_gamma_mu()
    test_q_epsilon()
    init_cotas()
