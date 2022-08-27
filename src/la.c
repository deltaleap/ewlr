#include <math.h>

#include "la.h"

void v2_copy(v2 *dest, v2 *src){
    /* copy values from dest vec (2d) to src vec (2d)
    */
   dest->x = src->x;
   dest->y = src->y;
}

void v3_copy(v3 *dest, v3 *src){
    /* copy values from dest vec (3d) to src vec (3d)
    */
   dest->x = src->x;
   dest->y = src->y;
   dest->z = src->z;
}

void v4_copy(v4 *dest, v4 *src){
    /* copy values from dest vec (4d) to src vec (4d)
    */
   dest->w = src->w;
   dest->x = src->x;
   dest->y = src->y;
   dest->z = src->z;
}

void m2_copy(m2 *dest, m2 *src){
    /* copy values from dest matrix (2x2) to src matrix (2x2)
    */
   dest->M[0][0] = src->M[0][0];
   dest->M[0][1] = src->M[0][1];
   dest->M[1][0] = src->M[1][0];
   dest->M[1][1] = src->M[1][1];
}

void m3_copy(m3 *dest, m3 *src){
    /* copy values from dest matrix (3x3) to src matrix (3x3)
    */
   for_loop_c(3)
       for_loop_r(3)
           dest->M[r][c] = src->M[r][c];
}

void m4_copy(m4 *dest, m4 *src){
    /* copy values from dest matrix (4x4) to src matrix (4x4)
    */
   for_loop_c(4)
       for_loop_r(4)
         dest->M[r][c] = src->M[r][c];

}

/* zero-vectors */
v2 v2_zero()
{
    /* TODO: test */
    v2 zero;
    zero.x = 0.f;
    zero.y = 0.f;
    return(zero);
}

v3 v3_zero()
{
    /* TODO: test */
    v3 zero;
    zero.x = 0.f;
    zero.y = 0.f;
    zero.z = 0.f;
    return(zero);
}

v4 v4_zero()
{
    /* TODO: test */
    v4 zero;
    zero.w = 0.f;
    zero.x = 0.f;
    zero.y = 0.f;
    zero.z = 0.f;
    return(zero);
}

/* zero-vectors: pointer version */
void v2_zero_(v2 *dest)
{
    /* TODO: test */
    dest->x = 0.f;
    dest->y = 0.f;
}

void v3_zero_(v3 *dest)
{
    /* TODO: test */
    dest->x = 0.f;
    dest->y = 0.f;
    dest->z = 0.f;
}

void v4_zero_(v4 *dest)
{
    /* TODO: test */
    dest->w = 0.f;
    dest->x = 0.f;
    dest->y = 0.f;
    dest->z = 0.f;
}

/* identity matrix */
m2 m2_eye()
{
    m2 eye;
    for_loop_r(2)
        for_loop_c(2)
            eye.M[r][c] = (r == c) ? 1.f : 0.f;
    return(eye);
}

m3 m3_eye()
{
    m3 eye;
    for_loop_r(3)
        for_loop_c(3)
            eye.M[r][c] = (r == c) ? 1.f : 0.f;
    return(eye);
}

m4 m4_eye()
{
    m4 eye;
    for_loop_r(4)
        for_loop_c(4)
            eye.M[r][c] = (r == c) ? 1.f : 0.f;
    return(eye);
}

/* identity matrix - pointer version */
void m2_eye_(m2 *dest)
{
    /* TODO: test */
    for_loop_r(2)
        for_loop_c(2)
            dest->M[r][c] = (r == c) ? 1.f : 0.f;
}

void m3_eye_(m3 *dest)
{
    /* TODO: test */
    for_loop_r(3)
        for_loop_c(3)
            dest->M[r][c] = (r == c) ? 1.f : 0.f;
}

void m4_eye_(m4 *dest)
{
    /* TODO: test */
    for_loop_r(4)
        for_loop_c(4)
            dest->M[r][c] = (r == c) ? 1.f : 0.f;
}

f32 v2_dot(v2 a, v2 b)
{
    /* dot product of 2d vectors
    */
    f32 res = a.x * b.x + a.y * b.y;
    return(res);
}

f32 v3_dot(v3 a, v3 b)
{
    /* dot product of 3d vectors
    */
    f32 res = a.x * b.x + a.y * b.y + a.z * b.z;
    return(res);
}

f32 v4_dot(v4 a, v4 b)
{
    /* dot product of 4d vectors
    */
    f32 res = a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
    return(res);
}

f32 m2_det(m2 mat)
{
    /* calculate determinant of 2x2 matrix
    */
    f32 det = mat.M[0][0] * mat.M[1][1] - mat.M[0][1] * mat.M[1][0];
    return(det);
}

f32 m3_det(m3 mat)
{
    /* calculate determinant of 3x3 matrix
    */
    m2 a = {{
        { mat.M[1][1], mat.M[1][2] },
        { mat.M[2][1], mat.M[2][2] }
    }};
    m2 b = {{
        { mat.M[1][0], mat.M[1][2] },
        { mat.M[2][0], mat.M[2][2] }
    }};
    m2 c = {{
        { mat.M[1][0], mat.M[1][1] },
        { mat.M[2][0], mat.M[2][1] }
    }};

    f32 det = mat.M[0][0] * m2_det(a) - mat.M[0][1] * m2_det(b) + mat.M[0][2] * m2_det(c);
    return(det);
}

f32 m4_det(m4 mat)
{
    /* calculate determinant of 4x4 matrix
    */
    m3 a = {{
        { mat.M[1][1], mat.M[1][2], mat.M[1][3] },
        { mat.M[2][1], mat.M[2][2], mat.M[2][3] },
        { mat.M[3][1], mat.M[3][2], mat.M[3][3] }
    }};

    m3 b = {{
        { mat.M[1][0], mat.M[1][2], mat.M[1][3] },
        { mat.M[2][0], mat.M[2][2], mat.M[2][3] },
        { mat.M[3][0], mat.M[3][2], mat.M[3][3] }
    }};

    m3 c = {{
        { mat.M[1][0], mat.M[1][1], mat.M[1][3] },
        { mat.M[2][0], mat.M[2][1], mat.M[2][3] },
        { mat.M[3][0], mat.M[3][1], mat.M[3][3] }
    }};

    m3 d = {{
        { mat.M[1][0], mat.M[1][1], mat.M[1][2] },
        { mat.M[2][0], mat.M[2][1], mat.M[2][2] },
        { mat.M[3][0], mat.M[3][1], mat.M[3][2] }
    }};

    f32 det = mat.M[0][0] * m3_det(a) - mat.M[0][1] * m3_det(b) + mat.M[0][2] * m3_det(c) - mat.M[0][3] * m3_det(d);
    return(det);
}

/* transpose */
m2 m2_transpose(m2 m)
{
    m2 res;
    for_loop_r(2)
        for_loop_c(2)
            res.M[r][c] = m.M[c][r];
    return(res);
}

m3 m3_transpose(m3 m)
{
    m3   res;
    for_loop_r(3)
        for_loop_c(3)
            res.M[r][c] = m.M[c][r];
    return(res);
}

m4 m4_transpose(m4 m)
{
    m4 res;
    for_loop_r(4)
        for_loop_c(4)
            res.M[r][c] = m.M[c][r];
    return(res);
}

/* cholesky decomposition */
m2 m2_llt(m2 m)
{
    m2 L;
    for_loop_r(2)
    {
        for_loop_c(2)
        {
            f32 sum = 0;
            for (int k = 0; k < c; ++k)
                sum += L.M[r][k] * L.M[c][k];

            if (r == c)
                L.M[r][c] = sqrt(m.M[r][r] - sum);
            else
            {
                if (r > c)
                    L.M[r][c] = (1.f / L.M[c][c] * (m.M[r][c] - sum));
                else
                    L.M[r][c] = 0;
            }
        }
    }
    return(L);
}

m3 m3_llt(m3 m)
{
    m3 L;
    for_loop_r(3)
    {
        for_loop_c(3)
        {
            f32 sum = 0;
            for (int k = 0; k < c; ++k)
                sum += L.M[r][k] * L.M[c][k];

            if (r == c)
                L.M[r][c] = sqrt(m.M[r][r] - sum);
            else
            {
                if (r > c)
                    L.M[r][c] = (1.f / L.M[c][c] * (m.M[r][c] - sum));
                else
                    L.M[r][c] = 0;
            }
        }
    }
    return(L);
}

m4 m4_llt(m4 m)
{
    m4 L;
    for_loop_r(4)
    {
        for_loop_c(4)
        {
            f32 sum = 0;
            for (int k = 0; k < c; ++k)
                sum += L.M[r][k] * L.M[c][k];

            if (r == c)
                L.M[r][c] = sqrt(m.M[r][r] - sum);
            else
            {
                if (r > c)
                    L.M[r][c] = (1.f / L.M[c][c] * (m.M[r][c] - sum));
                else
                    L.M[r][c] = 0;
            }
        }
    }
    return(L);
}

/* inverse of symmetric matrix */
m2 m2s_inv(m2 m)
{
    /* Inverse of A:
        - A is simmetric, semipositive definite
        - A A⁻¹ = I
        - A = LLᵀ
        - A⁻¹ = (LLᵀ)⁻¹ = L⁻ᵀL⁻¹
        - LLᵀA⁻¹ = I -> Lu = I, where u = LᵀA⁻¹
        - resolve Lu = I through forward substitution
        - resolve LᵀA⁻¹ = u through back substitution
    */
    /* cholesky decomposition -> L */
    m2 L = m2_llt(m);
    /* hard-coded forward substitution */
    m2 u;
    u.M[0][0] = 1.f / L.M[0][0];
    u.M[1][1] = 1.f / L.M[1][1];
    u.M[1][0] = - L.M[1][0] / ( L.M[0][0] * L.M[1][1] );
    /* hard-coded backward substitution */
    m2 res;
    res.M[1][0] = u.M[1][0] / L.M[1][1];
    res.M[1][1] = u.M[1][1] / L.M[1][1];
    /* not useful
    res.M[0][0] = ( u.M[0][0] - u.M[1][0] * L.M[1][0] / L.M[1][1] ) / L.M[0][0];
    res.M[0][1] = - ( u.M[1][1] * L.M[1][0] ) / ( L.M[0][0] * L.M[1][1] );
    */
    res.M[0][1] = res.M[1][0];
    res.M[0][0] = res.M[1][1];

    return(res);
}


m3 m3s_inv(m3 m)
{
    /* Inverse of A:
        - A is simmetric, semipositive definite
        - A A⁻¹ = I
        - A = LLᵀ
        - A⁻¹ = (LLᵀ)⁻¹ = L⁻ᵀL⁻¹
        - LLᵀA⁻¹ = I -> Lu = I, where u = LᵀA⁻¹
        - resolve Lu = I through forward substitution
        - resolve LᵀA⁻¹ = u through back substitution
    */
    m3 L = m3_llt(m);
    f32 u00, u10, u20, u11, u21, u22;
    u00 = 1.f / L.M[0][0];
    u11 = 1.f / L.M[1][1];
    u22 = 1.f / L.M[2][2];
    u10 = - u00 * L.M[1][0] / L.M[1][1];
    u21 = - u11 * L.M[2][1] / L.M[2][2];
    u20 = - (u00 * L.M[2][0] + u10 * L.M[2][1]) / L.M[2][2];

    f32 a00, a10, a20, a11, a21, a22;
    a20 = u20 / L.M[2][2];
    a21 = u21 / L.M[2][2];
    a22 = u22 / L.M[2][2];
    a11 = (u11 - a21 * L.M[2][1]) / L.M[1][1];
    a10 = (u10 - a20 * L.M[2][1]) / L.M[1][1];
    a00 = (u00 - a10 * L.M[1][0] - a20 * L.M[2][0]) / L.M[0][0];

    m3 res;
    res.M[0][0] = a00;
    res.M[1][0] = a10;
    res.M[2][0] = a20;
    res.M[1][1] = a11;
    res.M[2][1] = a21;
    res.M[2][2] = a22;
    res.M[0][1] = a10;
    res.M[0][2] = a20;
    res.M[1][2] = a21;

    return(res);
}

m4 m4s_inv(m4 m)
{
    /* Inverse of A:
        - A is simmetric, semipositive definite
        - A A⁻¹ = I
        - A = LLᵀ
        - A⁻¹ = (LLᵀ)⁻¹ = L⁻ᵀL⁻¹
        - LLᵀA⁻¹ = I -> Lu = I, where u = LᵀA⁻¹
        - resolve Lu = I through forward substitution
        - resolve LᵀA⁻¹ = u through back substitution
    */
    m4 L = m4_llt(m);
    f32 u00, u10, u11, u20, u21, u22, u30, u31, u32, u33;
    f32 a00, a10, a11, a20, a21, a22, a30, a31, a32, a33;
    u00 = 1.f / L.M[0][0];
    u11 = 1.f / L.M[1][1];
    u22 = 1.f / L.M[2][2];
    u33 = 1.f / L.M[3][3];
    u32 = - u22 * L.M[3][2] / L.M[3][3];
    u21 = - u11 * L.M[2][1] / L.M[2][2];
    u10 = - u00 * L.M[1][0] / L.M[1][1];
    u20 = - (u00 * L.M[2][0] + u10 * L.M[2][1]) / L.M[2][2];
    u31 = - (u11 * L.M[3][1] + u21 * L.M[3][2]) / L.M[3][3];
    u30 = - (u00 * L.M[3][0] + u10 * L.M[3][1] + u20 * L.M[3][2]) / L.M[3][3];
    a30 = u30 / L.M[3][3];
    a33 = u33 / L.M[3][3];
    a32 = u32 / L.M[3][3];
    a31 = u31 / L.M[3][3];
    a22 = ( u22 - L.M[3][2] * a32 ) / L.M[2][2];
    a21 = ( u21 - L.M[3][2] * a31 ) / L.M[2][2];
    a20 = ( u20 - L.M[3][2] * a30 ) / L.M[2][2];
    a11 = ( u11 - L.M[2][1] * a21 - L.M[3][1] * a31 ) / L.M[1][1];
    a10 = ( u10 - L.M[2][1] * a20 - L.M[3][1] * a30 ) / L.M[1][1];
    a00 = ( u00 - L.M[1][0] * a10 - L.M[2][0] * a20 - L.M[3][0] * a30 ) / L.M[0][0];

    m4 res;
    res.M[0][0] = a00;
    res.M[1][0] = a10;
    res.M[2][0] = a20;
    res.M[3][0] = a30;
    res.M[1][1] = a11;
    res.M[2][1] = a21;
    res.M[3][1] = a31;
    res.M[2][2] = a22;
    res.M[3][2] = a32;
    res.M[3][3] = a33;

    res.M[0][1] = a10;
    res.M[0][2] = a20;
    res.M[0][3] = a30;
    res.M[1][2] = a21;
    res.M[1][3] = a31;
    res.M[2][3] = a32;

    return(res);
}


/* scalar sum */
m2 m2_ssum(f32 a, m2 m)
{
    /* sum of a scalar to a 2x2 matrix
    */
    m2 res;
    res.M[0][0] = a + m.M[0][0];
    res.M[0][1] = a + m.M[0][1];
    res.M[1][0] = a + m.M[1][0];
    res.M[1][1] = a + m.M[1][1];
    return(res);
}

v2 v2_ssum(f32 a, v2 v)
{
    /* sum of a scalar to a 2d vector
    */
    v2 res;
    res.x = a + v.x;
    res.y = a + v.y;
    return(res);
}

m3 m3_ssum(f32 a, m3 m)
{
    /* sum of a scalar to a 3x3 matrix
    */
    m3 res;
    for (int c = 0; c < 3; ++c)
    {
        for (int r = 0; r < 3; ++r)
        {
            res.M[r][c] = a + m.M[r][c];
        }
    }
    return(res);
}

v3 v3_ssum(f32 a, v3 v)
{
    /* sum of a scalar to a 3d vector
    */
    v3 res;
    res.x = a + v.x;
    res.y = a + v.y;
    res.z = a + v.z;
    return(res);
}

m4 m4_ssum(f32 a, m4 m)
{
    /* sum of a scalar to a 3x3 matrix
    */
    m4 res;

    for (int c = 0; c < 4; ++c)
    {
        for (int r = 0; r < 4; ++r)
        {
            res.M[r][c] = a + m.M[r][c];
        }
    }
    return(res);
}

v4 v4_ssum(f32 a, v4 v)
{
    /* sum of a scalar to a 4d vector
    */
    v4 res;
    res.w = a + v.w;
    res.x = a + v.x;
    res.y = a + v.y;
    res.z = a + v.z;
    return(res);
}

/* matrix sum */
m2 m2_msum(m2 a, m2 b)
{
    /* sum between 2x2 matrices
    */
   m2 res;
   res.M[0][0] = a.M[0][0] + b.M[0][0];
   res.M[0][1] = a.M[0][1] + b.M[0][1];
   res.M[1][0] = a.M[1][0] + b.M[1][0];
   res.M[1][1] = a.M[1][1] + b.M[1][1];
   return(res);
}

v2 v2_msum(v2 a, v2 b)
{
    /* sum between 2d vectors
    */
   v2 res;
   res.x = a.x + b.x;
   res.y = a.y + b.y;
   return(res);
}

m3 m3_msum(m3 a, m3 b)
{
    /* sum between 3x3 matrices
    */
    m3 res;

    for (int c = 0; c < 3; ++c)
    {
        for (int r = 0; r < 3; ++r)
        {
            res.M[r][c] = a.M[r][c] + b.M[r][c];
        }
    }
    return(res);
}

v3 v3_msum(v3 a, v3 b)
{
    /* sum between 3d vectors
    */
   v3 res;
   res.x = a.x + b.x;
   res.y = a.y + b.y;
   res.z = a.z + b.z;
   return(res);
}

m4 m4_msum(m4 a, m4 b)
{
    /* sum between 4x4 matrices
    */
    m4 res;
    for (int c = 0; c < 4; ++c)
    {
        for (int r = 0; r < 4; ++r)
        {
            res.M[r][c] = a.M[r][c] + b.M[r][c];
        }
    }
    return(res);
}

v4 v4_msum(v4 a, v4 b)
{
    /* sum between 4d vectors
    */
   v4 res;
   res.w = a.w + b.w;
   res.x = a.x + b.x;
   res.y = a.y + b.y;
   res.z = a.z + b.z;
   return(res);
}

/* matrix substraction */
v2 v2_msub(v2 a, v2 b)
{
    v2 res;
    res.x = a.x - b.x;
    res.y = a.y - b.y;
    return(res);
}

v3 v3_msub(v3 a, v3 b)
{
    v3 res;
    res.x = a.x - b.x;
    res.y = a.y - b.y;
    res.z = a.z - b.z;
    return(res);
}

v4 v4_msub(v4 a, v4 b)
{
    v4 res;
    res.w = a.w - b.w;
    res.x = a.x - b.x;
    res.y = a.y - b.y;
    res.z = a.z - b.z;
    return(res);
}

m2 m2_msub(m2 a, m2 b)
{
    m2 res;
    start_matrix_loop(2)
        res.M[r][c] = a.M[r][c] - b.M[r][c];
    end_matrix_loop
    return(res);
}

m3 m3_msub(m3 a, m3 b)
{
    m3 res;
    start_matrix_loop(3)
        res.M[r][c] = a.M[r][c] - b.M[r][c];
    end_matrix_loop
    return(res);
}

m4 m4_msub(m4 a, m4 b)
{
    m4 res;
    start_matrix_loop(4)
        res.M[r][c] = a.M[r][c] - b.M[r][c];
    end_matrix_loop
    return(res);
}

/* scalar multiplication */
m2 m2_smul(f32 a, m2 m)
{
    /* multiplication between scalar and 2x2 matrices
    */
    m2 res;
    res.M[0][0] = a * m.M[0][0];
    res.M[0][1] = a * m.M[0][1];
    res.M[1][0] = a * m.M[1][0];
    res.M[1][1] = a * m.M[1][1];
    return(res);
}
v2 v2_smul(f32 a, v2 v)
{
    /* multiplication between scalar and 2d vector
    */
    v2 res;
    res.x = a * v.x;
    res.y = a * v.y;
    return(res);
}
m3 m3_smul(f32 a, m3 m)
{
    /* multiplication between scalar and 3x3 matrices
    */
    m3 res;

    for (int c = 0; c < 3; ++c)
    {
        for (int r = 0; r < 3; ++r)
        {
            res.M[r][c] = a * m.M[r][c];
        }
    }
    return(res);
}

v3 v3_smul(f32 a, v3 v)
{
    /* multiplication between scalar and 3d vector
    */
    v3 res;
    res.x = a * v.x;
    res.y = a * v.y;
    res.z = a * v.z;
    return(res);
}

m4 m4_smul(f32 a, m4 m)
{
    /* multiplication between scalar and 4x4 matrices
    */
    m4 res;

    for (int c = 0; c < 4; ++c)
    {
        for (int r = 0; r < 4; ++r)
        {
            res.M[r][c] = a * m.M[r][c];
        }
    }
    return(res);
}

v4 v4_smul(f32 a, v4 v)
{
    /* multiplication between scalar and 4d vector
    */
    v4 res;
    res.w = a * v.w;
    res.x = a * v.x;
    res.y = a * v.y;
    res.z = a * v.z;
    return(res);
}

/* scalar division */
v2 v2_sdiv(f32 a, v2 v)
{
    v2 res;
    res.x = v.x / a;
    res.y = v.y / a;
    return(res);
}

v3 v3_sdiv(f32 a, v3 v)
{
    v3 res;
    res.x = v.x / a;
    res.y = v.y / a;
    res.z = v.z / a;
    return(res);
}

v4 v4_sdiv(f32 a, v4 v)
{
    v4 res;
    res.w = v.w / a;
    res.x = v.x / a;
    res.y = v.y / a;
    res.z = v.z / a;
    return(res);
}

m2 m2_sdiv(f32 a, m2 m)
{
    m2 res;
    start_matrix_loop(2)
        res.M[r][c] = m.M[r][c] / a;
    end_matrix_loop
    return(res);
}

m3 m3_sdiv(f32 a, m3 m)
{
    m3 res;
    start_matrix_loop(3)
        res.M[r][c] = m.M[r][c] / a;
    end_matrix_loop
    return(res);
}
m4 m4_sdiv(f32 a, m4 m)
{
    m4 res;
    start_matrix_loop(4)
        res.M[r][c] = (m.M[r][c] / a);
    end_matrix_loop
    return(res);
}

/* matrix multiplication */
m2 m2_mmul(m2 a, m2 b)
{
    /* matrix multiplication between 2x2 matrices
        2x2, 2x2 -> 2x2 
    */
    m2 res;
    res.M[0][0] = a.M[0][0] * b.M[0][0] + a.M[0][1] * b.M[1][0];
    res.M[0][1] = a.M[0][0] * b.M[0][1] + a.M[0][1] * b.M[1][1];
    res.M[1][0] = a.M[1][0] * b.M[0][0] + a.M[1][1] * b.M[1][0];
    res.M[1][1] = a.M[1][0] * b.M[0][1] + a.M[1][1] * b.M[1][1];
    return(res);
}

m2 v2_mmul(v2 a, v2 b)
{
    /* matrix multiplication between 2d vectors
        2x1, 1x2 -> 2x2
    */
   m2 res;
   res.M[0][0] = a.x * b.x;
   res.M[0][1] = a.x * b.y;
   res.M[1][0] = a.y * b.x;
   res.M[1][1] = a.y * b.y;
   return(res);
}

v2 m2v2_mmul(m2 a, v2 b)
{
    /* matrix multiplication between 2x2 matrix and 2d vector
        2x2, 2x1 -> 2x1
    */
   v2 res;
   res.x = a.M[0][0] * b.x + a.M[0][1] * b.y;
   res.y = a.M[1][0] * b.x + a.M[1][1] * b.y;
   return(res);
}

v2 v2m2_mmul(v2 a, m2 b)
{
    /* matrix multiplication between 2d vector and 2x2 matrix
        1x2, 2x2 -> 1x2
    */
   v2 res;
   res.x = a.x * b.M[0][0] + a.y * b.M[1][0];
   res.y = a.x * b.M[0][1] + a.y * b.M[1][1];
   return(res);
}

m3 m3_mmul(m3 a, m3 b)
{
    /* matrix multiplication between 3x3 matrices
        3x3, 3x3 -> 3x3
    */
    m3 res;
    for (int Col = 0; Col < 3; ++Col)
    {
        for (int Row = 0; Row < 3; ++Row)
        {
            f32 sum = 0;
            int K;
            for (K = 0; K < 3; ++K)
            {
                sum += a.M[Row][K] * b.M[K][Col];
            }
            res.M[Row][Col] = sum;
        }
    }
    return(res);
}

m3 v3_mmul(v3 a, v3 b)
{
    /* matrix multiplication between 2d vectors
        3x1, 1x3 -> 3x3
    */
   m3 res;
   #include <stdio.h>
   res.M[0][0] = a.x * b.x;
   res.M[0][1] = a.x * b.y;
   res.M[0][2] = a.x * b.z;
   res.M[1][0] = a.y * b.x;
   res.M[1][1] = a.y * b.y;
   res.M[1][2] = a.y * b.z;
   res.M[2][0] = a.z * b.x;
   res.M[2][1] = a.z * b.y;
   res.M[2][2] = a.z * b.z;
   return(res);
}

v3 m3v3_mmul(m3 a, v3 b)
{
    /* matrix multiplication between 2x2 matrix and 2d vector
        3x3, 3x1 -> 3x1
    */
    v3 res;
    res.x = a.M[0][0]*b.x + a.M[0][1]*b.y + a.M[0][2]*b.z;
    res.y = a.M[1][0]*b.x + a.M[1][1]*b.y + a.M[1][2]*b.z;
    res.z = a.M[2][0]*b.x + a.M[2][1]*b.y + a.M[2][2]*b.z;
    return(res);
}

v3 v3m3_mmul(v3 a, m3 b)
{
    /* matrix multiplication between 2d vector and 2x2 matrix
        1x3, 3x3 -> 1x3
    */
   v3 res;
   res.x = a.x*b.M[0][0] + a.y*b.M[1][0] + a.z*b.M[2][0];
   res.y = a.x*b.M[0][1] + a.y*b.M[1][1] + a.z*b.M[2][1];
   res.z = a.x*b.M[0][2] + a.y*b.M[1][2] + a.z*b.M[2][2];
   return(res);
}

m4 m4_mmul(m4 a, m4 b)
{
    /* matrix multiplication between 4x4 matrices
        4x4, 4x4 -> 4x4
    */
    m4 res;
    for (int Col = 0; Col < 4; ++Col)
    {
        for (int Row = 0; Row < 4; ++Row)
        {
            f32 sum = 0;
            int K;
            for (K = 0; K < 4; ++K)
            {
                sum += a.M[Row][K] * b.M[K][Col];
            }
            res.M[Row][Col] = sum;
        }
    }

    return(res);
}

m4 v4_mmul(v4 a, v4 b)
{
    /* matrix multiplication between 2d vectors
        4x1, 1x4 -> 4x4
    */
   m4 res;
   res.M[0][0] = a.w * b.w;
   res.M[0][1] = a.w * b.x;
   res.M[0][2] = a.w * b.y;
   res.M[0][3] = a.w * b.z;
   res.M[1][0] = a.x * b.w;
   res.M[1][1] = a.x * b.x;
   res.M[1][2] = a.x * b.y;
   res.M[1][3] = a.x * b.z;
   res.M[2][0] = a.y * b.w;
   res.M[2][1] = a.y * b.x;
   res.M[2][2] = a.y * b.y;
   res.M[2][3] = a.y * b.z;
   res.M[3][0] = a.z * b.w;
   res.M[3][1] = a.z * b.x;
   res.M[3][2] = a.z * b.y;
   res.M[3][3] = a.z * b.z;
   return(res);
}

v4 m4v4_mmul(m4 a, v4 b)
{
    /* matrix multiplication between 2x2 matrix and 2d vector
        4x4, 4x1 -> 4x1
    */
    v4 res;
    res.w = a.M[0][0]*b.w + a.M[0][1]*b.x + a.M[0][2]*b.y + a.M[0][3]*b.z;
    res.x = a.M[1][0]*b.w + a.M[1][1]*b.x + a.M[1][2]*b.y + a.M[1][3]*b.z;
    res.y = a.M[2][0]*b.w + a.M[2][1]*b.x + a.M[2][2]*b.y + a.M[2][3]*b.z;
    res.z = a.M[3][0]*b.w + a.M[3][1]*b.x + a.M[3][2]*b.y + a.M[3][3]*b.z;
    return(res);
}

v4 v4m4_mmul(v4 a, m4 b)
{
    /* matrix multiplication between 2d vector and 2x2 matrix
        1x4, 4x4 -> 1x4
    */
   v4 res;
   res.w = a.w*b.M[0][0] + a.x*b.M[1][0] + a.y*b.M[2][0] + a.z*b.M[3][0];
   res.x = a.w*b.M[0][1] + a.x*b.M[1][1] + a.y*b.M[2][1] + a.z*b.M[3][1];
   res.y = a.w*b.M[0][2] + a.x*b.M[1][2] + a.y*b.M[2][2] + a.z*b.M[3][2];
   res.z = a.w*b.M[0][3] + a.x*b.M[1][3] + a.y*b.M[2][3] + a.z*b.M[3][3];
   return(res);
}