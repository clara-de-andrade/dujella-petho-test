#ifndef __BASE_H__
#define __BASE_H__

#include <stdint.h>
#include <mpfr.h>

/* ************************************************************************* */

void nint_dist(mpfr_t dist, mpfr_t x, mpfr_prec_t prec);
void qconv(mpz_t q, mpfr_t eta, uint64_t t, mpfr_prec_t prec);
void findroot(mpfr_t root, uint64_t k, mpfr_prec_t prec);

/* ************************************************************************* */

void prog_bar(uint64_t i, uint64_t i_max);

/* ************************************************************************* */

#endif
