#!/usr/bin/python3

"""
Seja k >= 2. Em adição, seja n -> F(k, n) a sequência definida pela relação de
recorrência F(k, n) = F(k, n - 1) + ... + F(k, n - k), com os valores iniciais
F(k, -(k - 2)) = ... = F(k, 0) = 0 e F(k, 1) = 1. Por sua vez, seja alpha_k a
raíz real do polinômio P(k, x) = x**k - x**(k - 1) - ... - x - 1 contida no
intervalo entre 2 * (1 - 2**-k) e 2. Por fim, seja g_k definido em termos de
alpha_k por g_k = (alpha_k - 1) / (2 + (k + 1) * (alpha - 2)). Em particular:

    |F(k, n) - g_k * alpha_k ** (n - 1)| < 1/2                              (1)

Com base no trabalho de Ana Paula de Araújo Chaves, em sua tese de doutorado,
Equações Diofantinas Envolvendo Potências de Termos de Sequências Recorrentes,
verificamos a seguir que não existem soluções inteiras (m, n, k, s) da equação

    F(k, m) ** s + F(k, m + 1) ** s = F(k, n)                               (2)

quando 3 <= m <= 1394 e 3 <= k <= min(m, log s). (Famosamente, esta equação
tem solução quando (m, n, k, s) = (m, 2m + 1, 2, 2).) Para tal, usamos alguns
lemas desenvolvidos na tese, relacionando m, n, k e s, e o seguinte teorema:

    Teorema. (Dujella-Petho) Suponha que M é um inteiro positivo. Seja p/q um
    convergente da expansão em fração continuada do número irracional gamma
    tal que q > 6 * M; e seja epsilon = ||mu * q|| - M * ||gamma * q||, onde
    mu é um número real, e ||x|| denota a distância entre x e o inteiro mais
    próximo de x. Se epsilon > 0, então não existe solução para a desigualdade

        0 < m * gamma - n + mu < A * B ** (-m)

    em inteiros positivos m, n com

        log(A * q / epsilon) / log(B) <= m < M.

Efetivamente, o programa presente computa F(k, m + 1) e aproximações racionais
de alpha_k e g_k para 0 <= m <= 1394 e 3 <= k <= log s, computa aproximações
racionais de

    gamma(k, m) = log(alpha_k) / log(F(k, m + 1))                           (3)

e

    mu(k, m) = log(alpha_k / g_k) / log(F(k, m + 1))                        (4)

para 3 <= m <= 1394 e 3 <= k <= min(m, log s), e enfim computa o denominador
q(k, m)_t do t-ésimo convergente de gamma(k, m), assim como uma aproximação de

    epsilon(k, m) = ||mu(k, m) * q(k, m)|| - M * ||gamma(k, m) * q(k, m)||, (5)

e verifica que q(k, m)_t > 6 * M e epsilon(k, m) > 0 para todo 3 <= m <= 1394
e 3 <= k <= min(m, log s), onde M e t são parâmetros pré-computados. Por fim,
após a verificação, computamos log(A * q(k, m)_t / epsilon(k, m)_t) / log(B)
para todo 3 <= m <= 1394 e 3 <= k <= min(m, log s), para determinar o mínimo
deste valor e encontrar novas cotas para n, s e k, permitindo assim verificar
o resultado principal em tempo finito.
"""

from mpmath import nstr

from tests import (
    init_cache,
    init_kbonacci,
    init_alpha_g,
    init_gamma_mu,
    test_q_epsilon,
    init_cotas,
)

"""
Passo 1. Inicialização de parâmetros

No que segue, todos os números decimais serão computados com 2 ** prec bits de
precisão; por exemplo, alpha[k] é uma aproximação racional do valor de alpha_k
que coincide nos 2 ** prec primeiros dígitos em base 2; e qualquer aproximação
com valor absoluto menor que 2 ** prec é arredondado para 0. Experimentando com
os valores de prec, pomos prec = 4096 para garantir que o convergente p_t/q_t
da aproximação de gamma(k, m) obtida seja igual ao convergente verdadeiro de
gamma(k, m), para nossa escolha do parâmetro t, e que a aproximação obtida de
epsilon(k, m)_t seja não nula. Por sua vez, pomos t >= t_min suficientemente
grande (i.e. t_min = 700) para que q(k, m)_t > 6 * M e epsilon(k, m)_t > 0.
Estes valores de t_min e prec são experimentais e podem ser ajustados.
"""

cache = init_cache(prec = 4096, m_max = 1394)

print(
    "Cotas iniciais:",
    f"  n <= {nstr(cache["n_max"], 3)}",
    f"  s <= {nstr(cache["s_max"], 3)}",
    f"  k <= {int(cache["k_max"])}",
    "",
    sep="\n"
)

"""
Passo 2. Pré-computação de F, alpha e g, gamma e mu

Por um resultado da tese, se existem soluções inteiras (m, n, k, s) da equação
(2) com 3 <= m <= 1394 e 3 <= k <= min(m, log s), então pela equação (1), temos

    |F(k, m + 1) ** s - g_k * alpha_k ** (n - 1)| < 2 * F(k, m) ** s.

Portanto, pondo Lambda = g_k * alpha_k ** (n - 1) / F(k, m + 1) ** s - 1, temos

    |Lambda| < 2 / (1.65 ** s)

e, pelo teorema de Matveev, obtemos que s <= 4.7e14 * (m ** 5) * (log(m) ** 4).

Usamos esta cota e as relações k <= log s e n <= m * s + 2 para determinar as
cotas de n, s e k, em função da cota m <= 1394, no passo anterior. Neste passo,
por sua vez, calculamos F(k, m + 1), alpha_k e g_k, para todo 3 <= m <= 1394 e
3 <= k <= min(m, log s), assim como os valores de gamma(k, m) e mu(k, m), dados
por (3) e (4), pois tomando log Lambda e dividindo por log F(k, m + 1), temos

    0 < n * gamma(k, m) - s - mu(k, m) < 2.01 * 1.65 ** (1 / 1394).         (6)

Isto é necessário para aplicar o Teorema de Dujella-Petho no próximo passo.
"""

cache.update(init_kbonacci(cache))
cache.update(init_alpha_g(cache))
cache.update(init_gamma_mu(cache))

"""
Passo 3. Verificação das condições do Teorema de Dujella-Petho

Para aplicar o Teorema de Dujella-Petho, é necessário verificar duas condições,
q(k, m)_t > 6 * M e epsilon(k, m)_t > 0, com t >= t_min, para 3 <= m <= 1394 e
3 <= k <= min(m, log s), onde M é a cota superior de n, pois pelo enunciado do
Teorema isto implica que não existem soluções da equação (6) para n inteiro e

    log(A * q(k, m)_t / epsilon(k, m)_t) / log(B) <= n < M,

onde A = 2.01 e B = 1.65 ** (1 / 1394). Simultaneamente, computamos o mínimo de

    min floor(log(A * q(k, m)_t / epsilon(k, m)_t) / log(B)),               (7)

com t >= t_min, para 3 <= m <= 1394 e 3 <= k <= min(m, log s).
"""

cache.update(test_q_epsilon(cache))

"""
Passo 4. Computação das cotas obtidas por Dujella-Petho

Satisfeitas as condições do Teorema de Dujella-Petho, calculamos o mínimo N_min
de (7), com t >= t_min, para 3 <= m <= 1394 e 3 <= k <= min(m, log s), Este
mínimo N_min é uma cota superior para n melhor que a cota superior inicialmente
calculada, e com ela podemos obter cotas melhores para s e k também, pois
s <= (n - 2) * (m - 1) e k <= log s, donde segue que a cota superior s_max de s
é igual a floor((N_min - 2) / 2), e por fim, no caso 3 <= m <= 1394, é
necessário verificar as soluções de (2) apenas nos intervalos

                  4 <= m <= 1394,
                  3 <= k <= floor(log s_max),
          ceil(e^3) <= s <= s_max,
    (m - 1) * s + 2 <= n <= m * s + 2.
"""

cache.update(init_cotas(cache))

print(
    "Cotas finais:",
    f"  n <= {cache["n_max"]}",
    f"  s <= {cache["s_max"]}",
    f"  k <= {cache["k_max"]}",
    "",
    sep="\n"
)
