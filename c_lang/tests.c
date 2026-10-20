#include <stdlib.h>
#include <assert.h>

#include "base.h"

#include "tests.h"

// TODO: pull values from cache
#define m_min 3
#define k_min 3

/* ************************************************************************* */

// TODO: init cache and set parameters

/* ************************************************************************* */

mpz_t** kbonacci;

void init_kbonacci(uint64_t k_max, uint64_t m_max)
{
    printf(
        "Pré-computando F(k, m + 1) para\n"
        "\t3 <= m <= %u\n"
        "\t3 <= k <= %u\n",
        m_max,
        k_max
    );

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

    printf("\nPRONTO!\n\n");
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
    printf(
        "Pré-computando alpha_k e g_k para\n"
        "\t3 <= k <= %u\n"
        "com precisão de %d bits\n",
        k_max,
        prec
    );

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

    printf("\nPRONTO!\n\n");
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

mpfr_t** gamma;
mpfr_t** mu;

void init_gamma_mu(uint64_t k_max, uint64_t m_max, mpfr_prec_t prec)
{
    printf(
        "Pré-computando gamma(k, m) e mu(k, m) para\n"
        "\t3 <= m <= %u\n"
        "\t3 <= k <= min(m, %u)\n"
        "com precisão de %d bits\n",
        m_max,
        k_max,
        prec
    );

    gamma = calloc(k_max + 1, sizeof(mpfr_t*));
    mu    = calloc(k_max + 1, sizeof(mpfr_t*));

    mpfr_t F, log_F;
    mpfr_inits2(prec, F, log_F, NULL);

    mpfr_t log_alpha[k_max + 1];
    mpfr_t log_g[k_max + 1];
    for (uint64_t k = k_min; k <= k_max; k++)
    {
        gamma[k] = calloc(m_max + 1, sizeof(mpfr_t));
        mu[k]    = calloc(m_max + 1, sizeof(mpfr_t));
        
        mpfr_inits2(prec, log_alpha[k], log_g[k], NULL);
        mpfr_log(log_alpha[k], alpha[k], MPFR_RNDN);
        mpfr_log(log_g[k], g[k], MPFR_RNDN);

    }

    for (uint64_t m = m_min; m <= m_max; m++)
    {
        for (uint64_t k = k_min; k <= (m < k_max ? m : k_max); k++)
        {
            mpfr_set_z(F, kbonacci[k][m + 1], MPFR_RNDN);
            mpfr_log(log_F, F, MPFR_RNDN);

            mpfr_init2(gamma[k][m], prec);
            mpfr_div(gamma[k][m], log_alpha[k], log_F, MPFR_RNDN);

            mpfr_init2(mu[k][m], prec);
            mpfr_sub(mu[k][m], log_alpha[k], log_g[k], MPFR_RNDN);
            mpfr_div(mu[k][m], mu[k][m], log_F, MPFR_RNDN);

        }
        
        prog_bar(m, m_max);
    }

    mpfr_clears(F, log_F, NULL);
    for (uint64_t k = k_min; k <= k_max; k++)
        mpfr_clears(log_alpha[k], log_g[k], NULL);
    
    printf("\nPRONTO!\n\n");
}

void clear_gamma_mu(uint64_t k_max, uint64_t m_max)
{
    for (uint64_t m = m_min; m <= m_max; m++)
    for (uint64_t k = k_min; k <= (m < k_max ? m : k_max); k++)
        mpfr_clears(gamma[k][m], mu[k][m], NULL);

    for (uint64_t k = k_min; k <= k_max; k++)
    {
        free(gamma[k]);
        free(mu[k]);
    }
    free(gamma);
    free(mu);
}

/* ************************************************************************* */

mpz_t** q;
mpfr_t** epsilon;

void test_q_epsilon(
    uint64_t k_max,
    uint64_t m_max,
    mpz_t n_max,
    uint64_t t_min,
    mpfr_prec_t prec)
{
    // TODO: print mpz_t in scientitic notation (base)
    gmp_printf(
        "Testando q(k, m)_t > 6 * n epsilon(k, m)_t > 0 para\n"
        "\t3 <= m <= %u\n"
        "\t3 <= k <= min(m, %u)\n"
        "\tn <= %Zd\n"
        "\tt >= %u\n"
        "com precisão de %d bits\n",
        m_max,
        k_max,
        n_max,
        t_min,
        prec
    );

    q       = calloc(k_max + 1, sizeof(mpz_t*));
    epsilon = calloc(k_max + 1, sizeof(mpfr_t*));

    mpz_t target;
    mpz_init(target);
    mpz_mul_ui(target, n_max, 6);

    mpfr_t t1, t2, q_fr, n_max_fr;
    mpfr_inits2(prec, t1, t2, q_fr, n_max_fr, NULL);

    for (uint64_t k = k_min; k <= k_max; k++)
    {
        q[k]       = calloc(m_max + 1, sizeof(mpz_t));
        epsilon[k] = calloc(m_max + 1, sizeof(mpfr_t));
    }

    for (uint64_t m = m_min; m <= m_max; m++)
    {
        for (uint64_t k = k_min; k <= (m < k_max ? m : k_max); k++)
        {
            mpz_init(q[k][m]);
            qconv(q[k][m], gamma[k][m], t_min, prec);
            assert(mpz_cmp(q[k][m], target) > 0);

            mpfr_set_z(q_fr, q[k][m], MPFR_RNDN);

            mpfr_mul(t1, q_fr, mu[k][m], MPFR_RNDN);
            nint_dist(t1, t1, prec);
            assert(mpfr_sgn(t1) != 0);

            mpfr_mul(t2, q_fr, gamma[k][m], MPFR_RNDN);
            nint_dist(t2, t2, prec);
            assert(mpfr_sgn(t2) != 0);

            mpfr_set_z(n_max_fr, n_max, MPFR_RNDN);
            mpfr_mul(t2, t2, n_max_fr, MPFR_RNDN);

            mpfr_init2(epsilon[k][m], prec);
            mpfr_sub(epsilon[k][m], t1, t2, MPFR_RNDN);
            assert(mpfr_sgn(epsilon[k][m]) > 0);
        }
        
        prog_bar(m, m_max);
    }

    printf("\nPRONTO!\n\n");
}

void clear_q_epsilon(uint64_t k_max, uint64_t m_max)
{

    for (uint64_t m = m_min; m <= m_max; m++)
    for (uint64_t k = k_min; k <= (m < k_max ? m : k_max); k++)
    {
        mpz_clear(q[k][m]);
        mpfr_clear(epsilon[k][m]);
    }

    for (uint64_t k = k_min; k <= k_max; k++)
    {
        free(q[k]);
        free(epsilon[k]);
    }
    free(q);
    free(epsilon);
}

/* ************************************************************************* */
