#ifndef __TESTS_H__

#include <stdint.h>
#include <stdio.h>

#include <gmp.h>
#include <mpfr.h>

/* ************************************************************************* */

extern uint64_t m_max;

extern uint64_t k_max;
extern mpz_t s_max, n_max;

void init_cache(mpfr_prec_t prec_init, uint64_t m_max_init);

/* ************************************************************************* */

extern mpz_t** F;

void init_F();
void clear_F();

/* ************************************************************************* */

extern mpfr_t* alpha;
extern mpfr_t* g;

void init_alpha_g();
void clear_alpha_g();

/* ************************************************************************* */

extern mpfr_t** gamma;
extern mpfr_t** mu;

void init_gamma_mu();
void clear_gamma_mu();

/* ************************************************************************* */

extern mpz_t** q;
extern mpfr_t** epsilon;

void test_q_epsilon();
void clear_q_epsilon();

/* ************************************************************************* */

void init_cotas();
void clear_cotas();

/* ************************************************************************* */

#endif
