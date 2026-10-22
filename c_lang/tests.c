#include <stdlib.h>
#include <assert.h>

#include "base.h"

#include "tests.h"

// TODO: pull values from cache
#define m_min 3
#define k_min 3

/* ************************************************************************* */

uint64_t m_max;

mpz_t s_max, n_max;
uint64_t k_max;

void init_cache(mpfr_prec_t prec_init, uint64_t m_max_init)
{
    printf(
        "Inicializando cachê com parâmetros:\n"
        "\tprec  = %ld\n"
        "\tm_max = %llu\n"
        "\n",
        prec_init,
        m_max_init
    );

    m_max = m_max_init;

    // Precisão dos cálculos com números decimais no mpfr, em bits.
    mpfr_set_default_prec(prec_init);

    mpfr_t temp_fr, s_max_fr, n_max_fr, k_max_fr;
    mpfr_inits(temp_fr, s_max_fr, n_max_fr, k_max_fr, NULL);
    
    // Cota superior de s. Dada por s <= 4.7e14 * (m ** 5) * (log(m) ** 4).
    mpz_init(s_max);
    {
        // temp_fr <- m_max
        mpfr_set_ui(temp_fr, m_max, MPFR_RNDN);
        // temp_fr <- log(temp_fr)
        mpfr_log(temp_fr, temp_fr, MPFR_RNDN);
        // temp_fr <- temp_fr ** 4
        mpfr_pow_ui(temp_fr, temp_fr, 4, MPFR_RNDN);

        // s_max_fr <- m_max
        mpfr_set_ui(s_max_fr, m_max, MPFR_RNDN);
        // s_max_fr <- s_max_fr ** 5
        mpfr_pow_ui(s_max_fr, s_max_fr, 5, MPFR_RNDN);
        // s_max_fr <- s_max_fr * temp_fr
        mpfr_mul(s_max_fr, s_max_fr, temp_fr, MPFR_RNDN);

        // s_max_fr <- s_max_fr * 47
        mpfr_mul_ui(s_max_fr, s_max_fr, 47, MPFR_RNDN);
        // temp_fr <- 10
        mpfr_set_ui(temp_fr, 10, MPFR_RNDN);
        // temp_fr <- temp_fr ** 13
        mpfr_pow_ui(temp_fr, temp_fr, 13, MPFR_RNDN);
        // s_max_fr <- s_max_fr * temp_fr
        mpfr_mul(s_max_fr, s_max_fr, temp_fr, MPFR_RNDN);
 
        // s_max <- floor(s_max_fr)
        mpfr_get_z(s_max, s_max_fr, MPFR_RNDD);
    }

    // Cota superior de n. Dada por n <= m * s + 2.
    mpz_init(n_max);
    {
        // n_max <- s_max * m_max
        // n_max <- n_max + 2
        mpz_mul_ui(n_max, s_max, m_max);
        mpz_add_ui(n_max, n_max, 2);

        // n_max_fr <- n_max
        mpfr_set_z(n_max_fr, n_max, MPFR_RNDN);
    }

    // Cota superior de k. Dada por k <= log(s).
    {
        // k_max_fr <- log(s_max_fr)
        mpfr_log(k_max_fr, s_max_fr, MPFR_RNDN);
        // k_max <- floor(k_max_fr)
        k_max = mpfr_get_ui(k_max_fr, MPFR_RNDD);
    }

    {
        mpfr_printf(
            "Cotas iniciais:\n"
            "\tn <= %.2Re\n"
            "\ts <= %.2Re\n"
            "\tk <= %llu\n"
            "\n",
            n_max_fr,
            s_max_fr,
            k_max
        );
    }

   mpfr_clears(temp_fr, s_max_fr, n_max_fr, k_max_fr, NULL);
}

/* ************************************************************************* */

mpz_t** F;

void init_F()
{
    printf(
        "Pré-computando F(k, m + 1) para\n"
        "\t3 <= m <= %llu\n"
        "\t3 <= k <= %llu\n",
        m_max,
        k_max
    );

    F = calloc(k_max + 1, sizeof(mpz_t*));

    for (uint64_t k = k_min; k <= k_max; k++)
    {
        F[k] = calloc(m_max + 2, sizeof(mpz_t));

        mpz_t prev[k];
        for (uint64_t i = 0; i < k; i++)
            if (i < k - 1)
                mpz_init_set_ui(prev[i], 0);
            else
                mpz_init_set_ui(prev[i], 1);

        mpz_init_set_ui(F[k][0], 0);
        mpz_init_set_ui(F[k][1], 1);

        for (uint64_t m = 2; m <= m_max + 1; m++)
        {
            mpz_init_set_ui(F[k][m], 0);
            for (uint64_t i = 0; i < k; i++)
                mpz_add(F[k][m], F[k][m], prev[i]);

            for (uint64_t i = 0; i < k; i++)
                if (i < k - 1)
                    mpz_set(prev[i], prev[i + 1]);
                else
                    mpz_set(prev[i], F[k][m]);
        }

        for (uint64_t i = 0; i < k; i++)
            mpz_clear(prev[i]);

        prog_bar(k, k_max);
    }

    printf("\nPRONTO!\n\n");
}

void clear_F()
{
    for (uint64_t k = k_min; k <= k_max; k++)
    {
        for (uint64_t m = 0; m <= m_max + 1; m++)
            mpz_clear(F[k][m]);
        free(F[k]);
    }
    free(F);
}

/* ************************************************************************* */

mpfr_t* alpha;
mpfr_t* g;

void init_alpha_g()
{
    mpfr_prec_t prec = mpfr_get_default_prec();

    printf(
        "Pré-computando alpha_k e g_k para\n"
        "\t3 <= k <= %llu\n"
        "com precisão de %ld bits\n",
        k_max,
        prec
    );

    alpha = calloc(k_max + 1, sizeof(mpfr_t));
    g     = calloc(k_max + 1, sizeof(mpfr_t));

    mpfr_t temp;
    mpfr_init(temp);

    mpfr_t p, dp, step, err;
    mpfr_inits(p, dp, step, err, NULL);
    // err <- 2 ** -prec
    mpfr_set_ui_2exp(err, 1, -prec, MPFR_RNDN);

    mpfr_t num, den;
    mpfr_inits(num, den, NULL);
 
    for (uint64_t k = k_min; k <= k_max; k++)
    {
        mpfr_init(alpha[k]);
        // alpha[k] <- 2
        mpfr_set_ui(alpha[k], 2, MPFR_RNDN);
        do
        {
            // p <- (alpha[k] ** k) * (alpha[k] - 2) + 1
            {
                // p <- alpha[k] ** k
                mpfr_pow_ui(p, alpha[k], k, MPFR_RNDN);
                // temp <- alpha[k] - 2
                mpfr_sub_ui(temp, alpha[k], 2, MPFR_RNDN);
                // p <- p * temp
                mpfr_mul(p, p, temp, MPFR_RNDN);
                // p <- p + 1
                mpfr_add_ui(p, p, 1, MPFR_RNDN);
            }

            // dp <- (alpha[k] ** (k - 1)) * ((k + 1) * alpha[k] - 2)
            {
                // dp <- alpha[k] ** (k - 1)
                mpfr_pow_ui(dp, alpha[k], k - 1, MPFR_RNDN);
                // temp <- (k + 1) * alpha[k]
                mpfr_mul_ui(temp, alpha[k], k + 1, MPFR_RNDN);
                // temp <- temp - 2 * k
                mpfr_sub_ui(temp, temp, 2 * k, MPFR_RNDN);
                // dp <- dp * temp
                mpfr_mul(dp, dp, temp, MPFR_RNDN);
            }

            // step <- p / dp
            mpfr_div(step, p, dp, MPFR_RNDN);
            // alpha[k] <- alpha[k] - step
            mpfr_sub(alpha[k], alpha[k], step, MPFR_RNDN);
        } while (mpfr_cmpabs(step, err) >= 0);

        mpfr_init(g[k]);
        // g[k] <- (alpha[k] - 1) / (2 + (k + 1) * (alpha[k] - 2))
        {

            // num <- alpha[k] - 1
            mpfr_sub_ui(num, alpha[k], 1, MPFR_RNDN);
            // den <- alpha[k] - 2
            mpfr_sub_ui(den, alpha[k], 2, MPFR_RNDN);
            // den <- (k + 1) * den
            mpfr_mul_ui(den, den, k + 1, MPFR_RNDN);
            // den <- 2 + den
            mpfr_add_ui(den, den, 2, MPFR_RNDN);
            // g[k] <- num / den
            mpfr_div(g[k], num, den, MPFR_RNDN);

        }

        prog_bar(k, k_max);
    }

    mpfr_clear(temp);
    mpfr_clears(p, dp, step, err, NULL);

    mpfr_clears(num, den, NULL);

    printf("\nPRONTO!\n\n");
}

void clear_alpha_g()
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

void init_gamma_mu()
{
    mpfr_prec_t prec = mpfr_get_default_prec();

    printf(
        "Pré-computando gamma(k, m) e mu(k, m) para\n"
        "\t3 <= m <= %llu\n"
        "\t3 <= k <= min(m, %llu)\n"
        "com precisão de %ld bits\n",
        m_max,
        k_max,
        prec
    );

    gamma = calloc(k_max + 1, sizeof(mpfr_t*));
    mu    = calloc(k_max + 1, sizeof(mpfr_t*));

    mpfr_t log_alpha[k_max + 1];
    mpfr_t log_g[k_max + 1];
    for (uint64_t k = k_min; k <= k_max; k++)
    {
        gamma[k] = calloc(m_max + 1, sizeof(mpfr_t));
        mu[k]    = calloc(m_max + 1, sizeof(mpfr_t));
        
        mpfr_inits(log_alpha[k], log_g[k], NULL);
        // log_alpha[k] <- log(alpha[k])
        mpfr_log(log_alpha[k], alpha[k], MPFR_RNDN);
        // log_g[k] <- log(g[k])
        mpfr_log(log_g[k], g[k], MPFR_RNDN);
    }

    mpfr_t F_fr, log_F;
    mpfr_inits(F_fr, log_F, NULL);

    for (uint64_t m = m_min; m <= m_max; m++)
    {
        for (uint64_t k = k_min; k <= (m < k_max ? m : k_max); k++)
        {
            // F_fr <- F(k, m + 1)
            mpfr_set_z(F_fr, F[k][m + 1], MPFR_RNDN);
            // log_F <- log(F_fr)
            mpfr_log(log_F, F_fr, MPFR_RNDN);

            mpfr_init(gamma[k][m]);
            // gamma[k][m] <- log_alpha[k] / log_F
            mpfr_div(gamma[k][m], log_alpha[k], log_F, MPFR_RNDN);

            mpfr_init(mu[k][m]);
            // mu[k][m] <- log_alpha[k] - log_g[k]
            mpfr_sub(mu[k][m], log_alpha[k], log_g[k], MPFR_RNDN);
            // mu[k][m] <- mu[k][m] / log_F
            mpfr_div(mu[k][m], mu[k][m], log_F, MPFR_RNDN);

        }
        
        prog_bar(m, m_max);
    }

    mpfr_clears(F_fr, log_F, NULL);

    for (uint64_t k = k_min; k <= k_max; k++)
        mpfr_clears(log_alpha[k], log_g[k], NULL);

    printf("\nPRONTO!\n\n");
}

void clear_gamma_mu()
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

void test_q_epsilon()
{
    mpfr_prec_t prec = mpfr_get_default_prec();

    mpfr_t n_max_fr;
    mpfr_init(n_max_fr);

    // n_max_fr <- n_max
    mpfr_set_z(n_max_fr, n_max, MPFR_RNDN);

    // TODO: print mpz_t in scientitic notation (base)
    mpfr_printf(
        "Testando q(k, m) > 6 * n epsilon(k, m) > 0 para\n"
        "\t3 <= m <= %llu\n"
        "\t3 <= k <= min(m, %llu)\n"
        "\tn <= %.2Re\n"
        "com precisão de %ld bits\n",
        m_max,
        k_max,
        n_max_fr,
        prec
    );

    q       = calloc(k_max + 1, sizeof(mpz_t*));
    epsilon = calloc(k_max + 1, sizeof(mpfr_t*));

    mpz_t target;
    mpz_init(target);
    // target <- n_max * 6
    mpz_mul_ui(target, n_max, 6);

    mpfr_t t1, t2, q_fr;
    mpfr_inits(t1, t2, q_fr, NULL);

    mpfr_t eta, a_fr;
    mpfr_inits(eta, a_fr, NULL);

    mpz_t q_old, temp, a;
    mpz_inits(q_old, temp, a, NULL);

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
            mpfr_init(epsilon[k][m]);

            // q_old <- 1
            mpz_set_ui(q_old, 1);
            // q[k][m] <- 0
            mpz_set_ui(q[k][m], 0);
            // eta <- gamma[k][m]
            mpfr_set(eta, gamma[k][m], MPFR_RNDN);

            while (1)
            {
                // a_fr <- floor(eta)
                mpfr_floor(a_fr, eta);
                // a <- a_fr
                mpfr_get_z(a, a_fr, MPFR_RNDN);

                // temp <- a * q[k][m]
                mpz_mul(temp, a, q[k][m]);
                // temp <- temp + q_old
                mpz_add(temp, temp, q_old);
                // q_old <- q[k][m]
                mpz_set(q_old, q[k][m]);
                // q[k][m] <- temp
                mpz_set(q[k][m], temp);

                // eta <- eta - a_fr
                mpfr_sub(eta, eta, a_fr, MPFR_RNDN);
                // eta <- 1 / eta
                mpfr_ui_div(eta, 1, eta, MPFR_RNDN);

                // if q[k][m] <= target
                if (mpz_cmp(q[k][m], target) <= 0)
                    continue;

                // q_fr <- q[k][m]
                mpfr_set_z(q_fr, q[k][m], MPFR_RNDN);

                // t1 <- q_fr * mu[k][m]
                mpfr_mul(t1, q_fr, mu[k][m], MPFR_RNDN);
                // t1 <- ||t1||
                nint_dist(t1, t1);

                // t2 <- q_fr * gamma[k][m]
                mpfr_mul(t2, q_fr, gamma[k][m], MPFR_RNDN);
                // t2 <- ||t2||
                nint_dist(t2, t2);


                // t2 <- t2 * n_max_fr
                mpfr_mul(t2, t2, n_max_fr, MPFR_RNDN);

                // epsilon[k][m] <- t1 - t2
                mpfr_sub(epsilon[k][m], t1, t2, MPFR_RNDN);

                // if epsilon[k][m] <= 0
                if (mpfr_sgn(epsilon[k][m]) <= 0)
                    continue;

                break;
            };
        }
        
        prog_bar(m, m_max);
    }

    mpz_clear(target);
    mpz_clears(q_old, temp, a, NULL);
    mpfr_clears(t1, t2, q_fr, n_max_fr, NULL);
    mpfr_clears(eta, a_fr, NULL);

    printf("\nPRONTO!\n\n");
}

void clear_q_epsilon()
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

void init_cotas()
{
    mpfr_prec_t prec = mpfr_get_default_prec();
    printf(
        "Reduzindo as cotas superiores de n, k, s\n"
        "com precisão de %ld bits\n",
        prec
    );

    mpfr_t log_A, log_B;
    mpfr_inits(log_A, log_B, NULL);

    {
        // log_A <- 2.01
        mpfr_set_str(log_A, "2.01", 10, MPFR_RNDN);
        // log_A <- log(log_A)
        mpfr_log(log_A, log_A, MPFR_RNDN);
    }

    {
        // log_B <- 1.65
        mpfr_set_str(log_B, "1.65", 10, MPFR_RNDN);
        // log_B <- log(log_B)
        mpfr_log(log_B, log_B, MPFR_RNDN);
    }

    mpfr_t temp_fr, n_max_fr, s_max_fr, k_max_fr, log_q, log_epsilon;
    mpfr_inits(temp_fr, n_max_fr, s_max_fr, k_max_fr, log_q, log_epsilon, NULL);

    // n_max_fr <- +inf
    mpfr_set_inf(n_max_fr, 1);

    for (uint64_t m = m_min; m <= m_max; m++)
    {
        for (uint64_t k = k_min; k <= (m < k_max ? m : k_max); k++)
        {
            // temp_fr <- q[k][m]
            mpfr_set_z(temp_fr, q[k][m], MPFR_RNDN);
            // log_q <- log(temp_fr)
            mpfr_log(log_q, temp_fr, MPFR_RNDN);
            // log_epsilon <- log(epsilon[k][m])
            mpfr_log(log_epsilon, epsilon[k][m], MPFR_RNDN);

            // temp_fr <- log_q - log_epsilon
            mpfr_sub(temp_fr, log_q, log_epsilon, MPFR_RNDN);
            // temp_fr <- temp_fr + log_A
            mpfr_add(temp_fr, temp_fr, log_A, MPFR_RNDN);
            // temp_fr <- temp_fr * m_max
            mpfr_mul_ui(temp_fr, temp_fr, m_max, MPFR_RNDN);
            // temp_fr <- temp_fr / log_B
            mpfr_div(temp_fr, temp_fr, log_B, MPFR_RNDN);

            // if temp_fr < n_max_fr
            if (mpfr_cmp(temp_fr, n_max_fr) < 0)
                // n_max_fr <- temp_fr
                mpfr_set(n_max_fr, temp_fr, MPFR_RNDN);
        }

        prog_bar(m, m_max);
    }

    // n_max <- n_max_fr
    mpfr_get_z(n_max, n_max_fr, MPFR_RNDD);

    {
        // s_max_fr <- n_max_fr - 2
        mpfr_sub_ui(s_max_fr, n_max_fr, 2, MPFR_RNDN);
        // s_max_fr <- s_max_fr / (m_min - 1)
        mpfr_div_ui(s_max_fr, s_max_fr, m_min - 1, MPFR_RNDN);
        // s_max <- s_max_fr
        mpfr_get_z(s_max, s_max_fr, MPFR_RNDD);
    }

    {
        // k_max_fr <- log(s_max_fr)
        mpfr_log(k_max_fr, s_max_fr, MPFR_RNDN);
        // k_max <- k_max_fr
        k_max = mpfr_get_ui(k_max_fr, MPFR_RNDD);
    }

    mpfr_clears(log_A, log_B, NULL);
    mpfr_clears(temp_fr, n_max_fr, s_max_fr,k_max_fr, log_q, log_epsilon, NULL);

    printf("\nPRONTO!\n\n");

    gmp_printf(
        "Cotas finais:\n"
        "\tn <= %Zd\n"
        "\ts <= %Zd\n"
        "\tk <= %llu\n",
        n_max,
        s_max,
        k_max
    );
}

void clear_cotas()
{
    mpz_clears(n_max, s_max, NULL);
}

/* ************************************************************************* */
