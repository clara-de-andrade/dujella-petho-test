#ifndef __BASE_H__
#define __BASE_H__

#include <stdint.h>
#include <mpfr.h>

/* ************************************************************************* */

typedef void (mpfr_fun_t)(mpfr_t, mpfr_t);

void findroot(mpfr_t root, uint64_t k, mpfr_t x0, mpfr_prec_t prec);

/* ************************************************************************* */

void prog_bar(uint64_t i, uint64_t i_max);

/* ************************************************************************* */

#endif