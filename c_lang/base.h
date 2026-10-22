#ifndef __BASE_H__
#define __BASE_H__

#include <stdint.h>
#include <mpfr.h>

/* ************************************************************************* */

void nint_dist(mpfr_t dist, mpfr_t x);

void prog_bar(uint64_t i, uint64_t i_max);

/* ************************************************************************* */

#endif
