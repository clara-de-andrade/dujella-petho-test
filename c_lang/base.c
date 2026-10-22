#include "base.h"

#include <stdio.h>

/* ************************************************************************* */

void nint_dist(mpfr_t dist, mpfr_t x)
{
    mpfr_t temp;
    mpfr_init(temp);
    
    mpfr_round(temp, x);
    mpfr_sub(temp, x, temp, MPFR_RNDN);
    mpfr_abs(dist, temp, MPFR_RNDN);

    mpfr_clear(temp);
}

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

