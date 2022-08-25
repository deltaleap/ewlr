#include <stdint.h>
#include <xmmintrin.h>

#define print2x2(M) printf("[[%f, %f]\n",M[0][0],M[0][1]); \
                    printf(" [%f, %f]]\n",M[1][0],M[1][1]);

#define print3x3(M) printf("[[%f, %f, %f]\n",M[0][0],M[0][1],M[0][2]); \
                    printf(" [%f, %f, %f]\n",M[1][0],M[1][1],M[1][2]); \
                    printf(" [%f, %f, %f]]\n",M[2][0],M[2][1],M[2][2]);

#define print4x4(M) printf("[[%f, %f, %f, %f]\n", M[0][0],M[0][1],M[0][2],M[0][3]); \
                    printf(" [%f, %f, %f, %f]\n", M[1][0],M[1][1],M[1][2],M[1][3]); \
                    printf(" [%f, %f, %f, %f]\n", M[2][0],M[2][1],M[2][2],M[2][3]); \
                    printf(" [%f, %f, %f, %f]]\n",M[3][0],M[3][1],M[3][2],M[3][3]);

#define for_loop_c(end) for(int c = 0; \
                                 c < end;  \
                                 ++c)
#define for_loop_r(end) for(int r = 0; \
                                 r < end;  \
                                 ++r)

#define start_matrix_loop(n) for(int c; c < n; ++c) \
                             { \
                                for(int r; r < n; ++r) \
                                {
#define end_matrix_loop }};

typedef float f32;
typedef uint8_t u8;

typedef struct v2 {
    f32 x, y;
} v2;
typedef struct v3 {
    f32 x, y, z;
} v3;
typedef struct v4 {
    f32 w, x, y, z;
} v4;
typedef struct m2 {
    f32 M[2][2];
} m2;
typedef struct m3 {
    f32 M[3][3];
} m3;
typedef struct m4 {
    f32 M[4][4];
} m4;

/* copy */
void v2_copy(v2 *dest, v2 *src);
void v3_copy(v3 *dest, v3 *src);
void v4_copy(v4 *dest, v4 *src);
void m2_copy(m2 *dest, m2 *src);
void m3_copy(m3 *dest, m3 *src);
void m4_copy(m4 *dest, m4 *src);

/* identity matrix */
m2 m2_eye();
m3 m3_eye();
m4 m4_eye();

/* dot product */
f32 v2_dot(v2 a, v2 b);
f32 v3_dot(v3 a, v3 b);
f32 v4_dot(v4 a, v4 b);

/* determinant */
f32 m2_det(m2 m);
f32 m3_det(m3 m);
f32 m4_det(m4 m);

/* transpose */
m2 m2_transpose(m2 m);
m3 m3_transpose(m3 m);
m4 m4_transpose(m4 m);

/* eigenvalues */
v2 m2_eigenval(m2 m);
v3 m3_eigenval(m3 m);
v4 m4_eigenval(m4 m);

/* eigenvectors */
m2 m2_eigenvec(m2 m);
m3 m3_eigenvec(m3 m);
m4 m4_eigenvec(m4 m);

/* cholesky decomposition */
m2 m2_llt(m2 m);
m3 m3_llt(m3 m);
m4 m4_llt(m4 m);

/* inverse of symmetric matrix */
m2 m2s_inv(m2 m);
m3 m3s_inv(m3 m);
m4 m4s_inv(m4 m);

/* inverse */
m2 m2_inv(m2 m);
m3 m3_inv(m3 m);
m4 m4_inv(m4 m);

/* scalar sum */
m2 m2_ssum(f32 a, m2 m);
v2 v2_ssum(f32 a, v2 v);
m3 m3_ssum(f32 a, m3 m);
v3 v3_ssum(f32 a, v3 v);
m4 m4_ssum(f32 a, m4 m);
v4 v4_ssum(f32 a, v4 v);

/* matrix sum */
m2 m2_msum(m2 a, m2 b);
v2 v2_msum(v2 a, v2 b);
m3 m3_msum(m3 a, m3 b);
v3 v3_msum(v3 a, v3 b);
m4 m4_msum(m4 a, m4 b);
v4 v4_msum(v4 a, v4 b);

/* scalar multiplication */
m2 m2_smul(f32 a, m2 m);
v2 v2_smul(f32 a, v2 v);
m3 m3_smul(f32 a, m3 m);
v3 v3_smul(f32 a, v3 v);
m4 m4_smul(f32 a, m4 m);
v4 v4_smul(f32 a, v4 v);

/* matrix multiplication */
m2 m2_mmul(m2 a, m2 b);   // 2x2, 2x2 -> 2x2 
m2 v2_mmul(v2 a, v2 b);   // 2x1, 1x2 -> 2x2
v2 m2v2_mmul(m2 a, v2 b); // 2x2, 2x1 -> 2x1
v2 v2m2_mmul(v2 a, m2 b); // 1x2, 2x2 -> 1x2

m3 m3_mmul(m3 a, m3 b);   // 3x3, 3x3 -> 3x3
m3 v3_mmul(v3 a, v3 b);   // 3x1, 1x3 -> 3x3
v3 m3v3_mmul(m3 a, v3 b); // 3x3, 3x1 -> 3x1
v3 v3m3_mmul(v3 a, m3 b); // 1x3, 3x3 -> 1x3

m4 m4_mmul(m4 a, m4 b);   // 4x4, 4x4 -> 4x4
m4 v4_mmul(v4 a, v4 b);   // 4x1, 1x4 -> 4x4
v4 m4v4_mmul(m4 a, v4 b); // 4x4, 4x1 -> 4x1
v4 v4m4_mmul(v4 a, m4 b); // 1x4, 4x4 -> 1x4