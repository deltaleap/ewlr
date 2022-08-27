/* Online Flexible least square
*/
#include "la.h"

typedef struct fls {
    f32 factor; // forgetting factor
    f32 mu; // (1 - factor) / factor
    v4 beta;
    m4 P, R;
} fls;

void fls_print  (fls *ls);
void fls_init   (fls *ls, f32 factor);
void fls_update (fls *ls, v4 x, f32 y);
f32  fls_predict(fls *ls, v4 x);