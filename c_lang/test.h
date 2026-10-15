#ifndef __TEST_H__
#define __TEST_H__

#include <stdint.h>
#include <stdio.h>

#include <gmp.h>
#include <mpfr.h>

/* ************************************************************************* */

extern mpz_t** kbonacci;
void init_kbonacci(uint64_t k_max, uint64_t m_max);
void clear_kbonacci(uint64_t k_max, uint64_t m_max);

/* ************************************************************************* */

extern mpfr_t* alpha;
extern mpfr_t* g;
void init_alpha_g(uint64_t k_max, mpfr_prec_t prec);
void clear_alpha_g(uint64_t k_max);

/* ************************************************************************* */

#endif