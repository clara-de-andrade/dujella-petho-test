#include "tests.h"
#include <locale.h>

#define k_max 77
#define m_max 1394
#define t_min 700

#define prec 4096

int main()
{
    #ifdef _WIN32
        setlocale(LC_ALL, ".UTF8");
    #else
        setlocale(LC_ALL, "");
    #endif

    init_kbonacci(k_max, m_max);
    init_alpha_g(k_max, prec);
    init_gamma_mu(k_max, m_max, prec);

    mpz_t n_max;
    mpz_init(n_max);
    mpz_set_ui(n_max, 10);
    mpz_pow_ui(n_max, n_max, 34);
    mpz_mul_ui(n_max, n_max, 941);
    test_q_epsilon(k_max, m_max, n_max, t_min, prec);

    clear_kbonacci(k_max, m_max);
    clear_alpha_g(k_max);
    clear_gamma_mu(k_max, m_max);
    clear_q_epsilon(k_max, m_max);

    return 0;
}
