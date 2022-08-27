/* Adaptive Forgetting EWRLS (AF-EWRLS)
*/

/* No new knowledge can be extracted from my telling.
   This confession has meant nothing.

   The AF-EWRLS algorithm (Algorithm 3, p. 75, Tsagaris 2010)
*/
#include "la.h"

typedef float f32;

typedef struct afewrls {
    v4 beta;
    v4 phi;
    m4 P;
    m4 S;
    f32 lambda;
    f32 lambda_minus;
    f32 lambda_plus;
    f32 alpha;
} afewrls;

void print_afewrls(afewrls *ls);
void afewrls_init(afewrls *ls, f32 alpha, f32 lambda_minus, f32 lambda_plus);
void afewrls_update(afewrls *ls, v4 x, f32 y);
f32 afewrls_predict(afewrls *ls, v4 x);