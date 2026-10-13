#include "base.h"

#include <stdio.h>

void findroot(mpfr_t root, uint64_t k, mpfr_t x0, mpfr_prec_t prec)
{

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