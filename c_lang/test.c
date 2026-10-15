#include <stdlib.h>

#include "base.h"

#include "test.h"

// TODO: pull values from cache
#define m_min 3
#define k_min 3

/* ************************************************************************* */

/* ************************************************************************* */

mpz_t** kbonacci;

void init_kbonacci(uint64_t k_max, uint64_t m_max)
{
    kbonacci = calloc(k_max + 1, sizeof(mpz_t*));

    for (uint64_t k = k_min; k <= k_max; k++)
    {
        kbonacci[k] = calloc(m_max + 2, sizeof(mpz_t));

        mpz_init_set_ui(kbonacci[k][0], 0);
        mpz_init_set_ui(kbonacci[k][1], 1);

        mpz_t prev[k];
        for (uint64_t i = 0; i < k; i++)
            if (i < k - 1)
                mpz_init_set_ui(prev[i], 0);
            else
                mpz_init_set_ui(prev[i], 1);

        mpz_t temp;
        mpz_init(temp);
        for (uint64_t m = 2; m <= m_max + 1; m++)
        {
            mpz_set_ui(temp, 0);
            for (uint64_t i = 0; i < k; i++)
                mpz_add(temp, temp, prev[i]);

            mpz_set(kbonacci[k][m], temp);

            for (uint64_t i = 0; i < k; i++)
                if (i < k - 1)
                    mpz_set(prev[i], prev[i + 1]);
                else
                    mpz_set(prev[i], temp);
        }

        for (uint64_t i = 0; i < k; i++)
            mpz_clear(prev[i]);
        mpz_clear(temp);

        prog_bar(k, k_max);
    }

    printf("\nPRONTO!\n");
}

void clear_kbonacci(uint64_t k_max, uint64_t m_max)
{
    for (uint64_t k = k_min; k <= k_max; k++)
    {
        for (uint64_t m = 0; m <= m_max + 1; m++)
            mpz_clear(kbonacci[k][m]);
        free(kbonacci[k]);
    }
    free(kbonacci);
}

/* ************************************************************************* */

mpfr_t* alpha;
mpfr_t* g;

void init_alpha_g(uint64_t k_max, mpfr_prec_t prec)
{
    alpha = calloc(k_max + 1, sizeof(mpfr_t));
    g     = calloc(k_max + 1, sizeof(mpfr_t));

    mpfr_t num, den;
    mpfr_inits2(prec, num, den, NULL);

    for (uint64_t k = k_min; k <= k_max; k++)
    {
        mpfr_init2(alpha[k], prec);
        findroot(alpha[k], k, prec);
        
        mpfr_sub_ui(num, alpha[k], 1, MPFR_RNDN);
        mpfr_sub_ui(den, alpha[k], 2, MPFR_RNDN);
        mpfr_mul_ui(den, den, k + 1, MPFR_RNDN);
        mpfr_add_ui(den, den, 2, MPFR_RNDN);

        mpfr_init2(g[k], prec);
        mpfr_div(g[k], num, den, MPFR_RNDN);

        prog_bar(k, k_max);
    }

    mpfr_clears(num, den, NULL);

    printf("\nPRONTO!\n");
}

void clear_alpha_g(uint64_t k_max)
{
    for (uint64_t k = k_min; k <= k_max; k++)
    {
        mpfr_clear(alpha[k]);
        mpfr_clear(g[k]);
    }
    free(alpha);
    free(g);
}

/* ************************************************************************* */