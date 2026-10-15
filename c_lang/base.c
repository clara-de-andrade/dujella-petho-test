#include "base.h"

#include <stdio.h>

/* ************************************************************************* */

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