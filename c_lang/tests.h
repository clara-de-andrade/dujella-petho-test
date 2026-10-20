#ifndef __TESTS_H__

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

extern mpfr_t** gamma;
extern mpfr_t** mu;
void init_gamma_mu(uint64_t k_max, uint64_t m_max, mpfr_prec_t prec);
void clear_gamma_mu(uint64_t k_max, uint64_t m_max);

/* ************************************************************************* */

extern mpz_t** q;
extern mpfr_t** epsilon;
extern mpfr_t** N;
void test_q_epsilon(
        uint64_t k_max,
        uint64_t m_max,
        mpz_t n_max,
        uint64_t t_min,
        mpfr_prec_t prec);
void clear_q_epsilon(uint64_t k_max, uint64_t m_max);

/* ************************************************************************* */

#endif
