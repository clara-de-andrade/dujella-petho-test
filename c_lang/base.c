#include "base.h"

#include <stdio.h>

/* ************************************************************************* */

void nint_dist(mpfr_t dist, mpfr_t x, mpfr_prec_t prec)
{
    mpfr_t temp;
    mpfr_init2(temp, prec);
    
    mpfr_round(temp, x);
    mpfr_sub(temp, x, temp, MPFR_RNDN);
    mpfr_abs(dist, temp, MPFR_RNDN);

    mpfr_clear(temp);
}

void qconv(mpz_t q, mpfr_t eta, uint64_t t_min, mpfr_prec_t prec)
{
    mpz_t u, v, temp, a;
    mpz_inits(u, v, temp, a, NULL);

    mpz_set_ui(u, 1);
    mpz_set_ui(v, 0);

    mpfr_t eta_new, a_fr;
    mpfr_inits2(prec, eta_new, a_fr, NULL);

    mpfr_set(eta_new, eta, MPFR_RNDN);

    for (uint64_t t = 0; t < t_min; t++)
    {
        mpfr_floor(a_fr, eta_new);
        mpfr_get_z(a, a_fr, MPFR_RNDN);

        mpz_mul(temp, a, v);
        mpz_add(temp, temp, u);
        mpz_set(u, v);
        mpz_set(v, temp);

        mpfr_sub(eta_new, eta_new, a_fr, MPFR_RNDN);
        mpfr_ui_div(eta_new, 1, eta_new, MPFR_RNDN);
    }

    mpz_set(q, v);

    mpz_clears(u, v, temp, a, NULL);
    mpfr_clears(eta_new, a_fr, NULL);
}

void findroot(mpfr_t root, uint64_t k, mpfr_prec_t prec)
{    
    mpfr_t x, step, p, dp, temp;
    mpfr_inits2(prec, x, step, p, dp, temp, NULL);

    {
        mpfr_set_ui_2exp(x, 1, 1, MPFR_RNDN);
        mpfr_set_ui_2exp(temp, 1, 1L - (long) k, MPFR_RNDN);
        mpfr_sub(x, x, temp, MPFR_RNDN);
        mpfr_set_ui_2exp(temp, k, -2L * (long) k, MPFR_RNDN);
        mpfr_mul_ui(temp, temp, k, MPFR_RNDN);
        mpfr_sub(x, x, temp, MPFR_RNDN);
    }

    do {
        mpfr_pow_ui(p, x, 2, MPFR_RNDN);
        mpfr_mul_2ui(temp, x, 1, MPFR_RNDN);
        mpfr_sub(p, p, temp, MPFR_RNDN);
        mpfr_pow_ui(temp, x, k - 1, MPFR_RNDN);
        mpfr_ui_div(temp, 1, temp, MPFR_RNDN);
        mpfr_add(p, p, temp, MPFR_RNDN);
        
        mpfr_mul_ui(dp, x, k + 1, MPFR_RNDN);
        mpfr_sub_ui(dp, dp, 2 * k, MPFR_RNDN);

        mpfr_div(step, p, dp, MPFR_RNDN);
        mpfr_sub(x, x, step, MPFR_RNDN);

        mpfr_abs(temp, step, MPFR_RNDN);
    } while (mpfr_cmp_ui_2exp(temp, 1, -(long)prec) >= 0);
    
    mpfr_set(root, x, MPFR_RNDN);

    mpfr_clears(x, step, p, dp, temp, NULL);
}

/* ************************************************************************* */

void prog_bar(uint64_t i, uint64_t i_max)
{
    int prog = 100 * i / i_max;
    int completed = prog / 2;
    int remaining = (100 - prog + 1) / 2;

    putchar('[');
    for (int i = 0; i < completed; i++)
        putchar('#');
    for (int i = 0; i < remaining; i++)
        putchar(' ');
    putchar(']');

    printf(" %i%%   \r", prog);
}
