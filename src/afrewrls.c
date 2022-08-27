#include "afrewrls.h"

#include <stdio.h>

void print_afewrls(afewrls *ls)
{
    printf("\nAFEWRLS[α=%.3f]\n", ls->alpha);
    printf("   β:\n   ");
    print4(ls->beta);
}

void afewrls_init(afewrls *ls, f32 alpha, f32 lambda_minus, f32 lambda_plus)
{
    v4_zero_(&(ls->beta));
    v4_zero_(&(ls->phi));
    m4_eye_(&(ls->P));
    m4_eye_(&(ls->S));
    ls->alpha = alpha;
    ls->lambda_minus = lambda_minus;
    ls->lambda_plus = lambda_plus;
    ls->lambda = 0.f;
}

void afewrls_update(afewrls *ls, v4 x, f32 y)
{
    v4 last_beta, last_phi;
    m4 last_S, last_P;
    v4_copy(&last_beta, &(ls->beta));
    v4_copy(&last_phi, &(ls->phi));
    m4_copy(&last_S, &(ls->S));
    m4_copy(&last_P, &(ls->P));

    m4 I = m4_eye();

    // calculte e
    f32 e = y - v4_dot(last_beta, x);
    printf("e: %.9f\n", e);

    // update lambda(t)
    f32 new_lambda;
    new_lambda = ls->lambda + ls->alpha * v4_dot(last_phi, x) * e;
    // truncate lambda(t)
    if (new_lambda > ls->lambda_plus)
        new_lambda = ls->lambda_plus;
    else if (new_lambda < ls->lambda_minus)
        new_lambda = ls->lambda_minus;
    ls->lambda = new_lambda;

    printf("new_lambda: %.9f\n", new_lambda);

    // update k(t)
    v4 k = v4_sdiv(
        (1 - v4_dot(x, m4v4_mmul(last_P, x)) / new_lambda),
        v4_sdiv(new_lambda, m4v4_mmul(last_P, x))
    );
    printf("k:\n");
    print4(k);

    // update beta(t)
    ls->beta = v4_msum(last_beta, v4_smul(e, k));
    printf("new_beta\n");
    print4(ls->beta);

    printf("new_lambda: %.9f\n", new_lambda);
    // update P(t)
    ls->P = m4_msub(
        m4_sdiv(new_lambda, last_P),
        m4_sdiv(new_lambda, m4_mmul(v4_mmul(k, x), last_P))
    );
    printf("P\n");
    print4x4(ls->P.M);

    printf("updating\n");
    // update S(t)
    m4 I_minus_kx = m4_msub(
        I,
        v4_mmul(k, x)
    );
    ls->S = m4_sdiv(
        new_lambda,
        m4_msub(
            m4_mmul(m4_mmul(I_minus_kx, last_S), I_minus_kx),
            m4_msum(ls->P, v4_mmul(k, k))
        )
    );
    printf("S\n");
    print4x4(ls->S.M);

    // update phi(t)
    ls->phi = v4_msum(
        m4v4_mmul(I_minus_kx, last_phi),
        v4_smul(e, m4v4_mmul(ls->S, x))
    );
    /* TODO: check correctness of phi:
       C version produces different results from py one.
       The error is in the second part of the formula (after the -)

       C: v4_smul(e, m4v4_mmul(ls->S, x))
       py: np.dot(new_S, x)
    */
    printf("phi part 1:\n");
    print4(m4v4_mmul(I_minus_kx, last_phi));
    printf("phi part 2:\n");
    print4(v4_smul(e, m4v4_mmul(ls->S, x)));
    printf("phi part 2.1:\n");
    print4(m4v4_mmul(ls->S, x));
    printf("phi\n");
    print4(ls->phi);
}

f32 afewrls_predict(afewrls *ls, v4 x)
{
   f32 res = v4_dot(ls->beta, x);
   return(res); 
}