/* Online Flexible least square
*/
#include <stdio.h>
#include "fls.h"

void fls_init(fls *ls, f32 factor)
{
    ls->factor = factor;
    v4_zero_(&(ls->beta));
    m4_eye_(&(ls->P));
    ls->mu = (1 - factor) / factor;
    ls->R = m4_sdiv(ls->mu, ls->P);
}

void fls_update (fls *ls, v4 x, f32 y)
{
    // compute vector k (with old R and x)
    v4 k = v4_sdiv(
        (1 + v4_dot(x, m4v4_mmul(ls->R, x))),
        m4v4_mmul(ls->R, x)
    );
    printf("k: \n");
    print4(k);
    // store old P (for updating R)
    m4 oldP;
    m4_copy(&(oldP), &(ls->P));
    // update P (with old R, k and x)
    ls->P = m4_msub(
        ls->R,
        m4_mmul(v4_mmul(k, x), ls->R)
    );
    // update R (with old P)
    ls->R = m4_sdiv(ls->mu, oldP);
    f32 e = y - v4_dot(ls->beta, x);
    ls->beta = v4_msum(ls->beta, v4_smul(e, k));
    print4(ls->beta);
    print4x4(ls->P.M);
    print4x4(ls->R.M);
}

f32  fls_predict(fls *ls, v4 x)
{
    return v4_dot(ls->beta, x);
}