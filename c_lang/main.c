#include "tests.h"
#include <locale.h>

#define prec  4096
#define m_max 1394

int main()
{
    #ifdef _WIN32
        setlocale(LC_ALL, ".UTF8");
    #else
        setlocale(LC_ALL, "");
    #endif

    init_cache(prec, m_max);

    init_F();
    init_alpha_g();
    init_gamma_mu();
    test_q_epsilon();
    init_cotas();

    clear_F();
    clear_alpha_g();
    clear_gamma_mu();
    clear_q_epsilon();
    clear_cotas();

    return 0;
}
