#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "la.h"
#include "afrewrls.h"
#include "fls.h"

#define EPSILON 0.0001

void test_la()
{
    printf("\n\ntesting dlMM linear algebra library.\n");
    printf("\n+ test **_copy +\n");

    printf(" - test v2\n");
    v2 s2 = { 1.0f, 1.2f };
    v2 d2;
    v2_copy(&d2, &s2);
    assert(s2.x == d2.x);
    assert(s2.y == d2.y);
    printf("     -> OK!\n");

    printf(" - test v3\n");
    v3 s3 = { 1.0f, 1.2f, 4.0f };
    v3 d3;
    v3_copy(&d3, &s3);
    assert(s3.x == d3.x);
    assert(s3.y == d3.y);
    assert(s3.z == d3.z);
    printf("     -> OK!\n");

    printf(" - test v4\n");
    v4 s4 = { 1.0f, 1.2f, 6.0f, 44.f };
    v4 d4;
    v4_copy(&d4, &s4);
    assert(s4.w == d4.w);
    assert(s4.x == d4.x);
    assert(s4.y == d4.y);
    assert(s4.z == d4.z);
    printf("     -> OK!\n");

    printf(" - test m2\n");
    m2 s22 = {{
        { 1.f, 2.f },
        { 3.f, 4.f }
    }};
    m2 d22;
    m2_copy(&d22, &s22);
    for (int c = 0; c < 2; ++c)
    {
        for (int r = 0; r < 2; ++r)
        {
            assert(s22.M[r][c] == d22.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf(" - test m3\n");
    m3 s33 = {{
        { 1.f, 2.f, 10.f },
        { 1.f, 2.f, 10.f },
        { 3.f, 4.f, 10.f }
    }};
    m3 d33;
    m3_copy(&d33, &s33);
    for (int c = 0; c < 3; ++c)
    {
        for (int r = 0; r < 3; ++r)
        {
            assert(s33.M[r][c] == d33.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf(" - test m4\n");
    m4 s44 = {{
        { 1.f, 2.f, 3.f, 4.f },
        { 11.f, 22.f, 33.f, 44.f },
        { 111.f, 222.f, 333.f, 444.f },
        { 1.23f, 2.34f, 3.45f, 4.56f }
    }};
    m4 d44;
    m4_copy(&d44, &s44);
    for (int c = 0; c < 4; ++c)
    {
        for (int r = 0; r < 4; ++r)
        {
            assert(s44.M[r][c] == d44.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf("\n+ test m*_eye +\n");

    printf(" - test m2\n");
    m2 eye2 = m2_eye();
    m2 expected_eye2 = {{
        { 1.f, 0.f },
        { 0.f, 1.f }
    }};
    for_loop_r(2)
        for_loop_c(2)
            assert(eye2.M[r][c] == expected_eye2.M[r][c]);
    printf("     -> OK!\n");

    printf(" - test m3\n");
    m3 eye3 = m3_eye();
    m3 expected_eye3 = {{
        { 1.f, 0.f, 0.f },
        { 0.f, 1.f, 0.f },
        { 0.f, 0.f, 1.f }
    }};
    for_loop_r(3)
        for_loop_c(3)
            assert(eye3.M[r][c] == expected_eye3.M[r][c]);
    printf("     -> OK!\n");

    printf(" - test m4\n");
    m4 eye4 = m4_eye();
    m4 expected_eye4 = {{
        { 1.f, 0.f, 0.f, 0.f },
        { 0.f, 1.f, 0.f, 0.f },
        { 0.f, 0.f, 1.f, 0.f },
        { 0.f, 0.f, 0.f, 1.f }
    }};
    for_loop_r(4)
        for_loop_c(4)
            assert(eye4.M[r][c] == expected_eye4.M[r][c]);
    printf("     -> OK!\n");

    printf("\n+ test v*_dot +\n");

    printf(" - test v2\n");
    v2 ar2 = { 10.f, 10.f };
    f32 rv2 = v2_dot(ar2, ar2);
    assert( rv2 == 200.f );
    printf("     -> OK!\n");

    printf(" - test v3\n");
    v3 ar3 = { 10.f, 10.f, 10.f };
    f32 rv3 = v3_dot(ar3, ar3);
    assert( rv3 == 300.f );
    printf("     -> OK!\n");

    printf(" - test v4\n");
    v4 ar4 = { 10.f, 10.f, 10.f, 10.f };
    f32 rv4 = v4_dot(ar4, ar4);
    assert( rv4 == 400.f );
    printf("     -> OK!\n");

    printf("\n+ test m*_det +\n");
    printf(" - test m2 (det = 0)\n");
    m2 mz2 = {{
        { 10.f, 10.f },
        { 10.f, 10.f }
    }};
    assert( m2_det(mz2) == 0.f );
    printf("     -> OK!\n");
    printf(" - test m2 (det /= 0)\n");
    m2 mnz2 = {{
        { 1.f, 2.f },
        { 3.f, 4.f }
    }};
    assert( m2_det(mnz2) == -2.f );
    printf("     -> OK!\n");

    printf(" - test m3\n");
    printf(" - test m3 (det = 0)\n");
    m3 mz3 = {{
        { 10.f, 10.f, 2.f },
        { 10.f, 10.f, 2000.f },
        { 10.f, 10.f, 2000.f }
    }};
    assert( m3_det(mz3) == 0.f );
    printf("     -> OK!\n");
    printf(" - test m3 (det /= 0)\n");
    m3 mnz3 = {{
        { 10.f, 210.f, 2.f },
        { 0.2f, 33.f, 12.f },
        { 2.1f, 99.f, 0.8f }
    }};
    assert( m3_det(mnz3) == -6456.6f );
    printf("     -> OK!\n");

    printf(" - test m4\n");
    printf(" - test m4 (det = 0)\n");
    m4 mz4 = {{
        { 10.f, 10.f, 2.f, 100.f },
        { 10.f, 10.f, 2000.f, 100.f  },
        { 10.f, 10.f, 2000.f, 100.f },
        { 10.f, 10.f, 2000.f, 100.f }
    }};
    assert( m4_det(mz4) == 0.f );
    printf("     -> OK!\n");
    printf(" - test m4 (det /= 0)\n");
    m4 mnz4 = {{
        { 10.f, 210.f, 2.f, 9.f },
        { 0.2f, 33.f, 12.f, 99.f },
        { 2.1f, 99.f, 0.8f, 1.2f },
        { 1.1f, 1.2f, 5.9f, 41.f }
    }};
    assert( m4_det(mnz4) - 60983.622f < 1.f );
    printf("     -> OK!\n");

    printf("\n+ test m*_transpose +\n");
    printf(" - test m2\n");
    m2 mtr2 = {{
        { 1.f, 3.f },
        { 2.f, 4.f }
    }};
    m2 rtr2 = m2_transpose(mtr2);
    m2 etr2 = {{
        { 1.f, 2.f },
        { 3.f, 4.f }
    }};
    for_loop_r(2)
        for_loop_c(2)
            assert(rtr2.M[r][c] == etr2.M[r][c]);
    printf("     -> OK!\n");

    printf(" - test m3\n");
    m3 mtr3 = {{
        { 1.f, 3.f, 5.f },
        { 2.f, 4.f, 6.f },
        { 2.f, 4.f, 6.f }
    }};
    m3 rtr3 = m3_transpose(mtr3);
    m3 etr3 = {{
        { 1.f, 2.f, 2.f },
        { 3.f, 4.f, 4.f },
        { 5.f, 6.f, 6.f }
    }};
    for_loop_r(3)
        for_loop_c(3)
            assert(rtr3.M[r][c] == etr3.M[r][c]);
    printf("     -> OK!\n");

    printf(" - test m4\n");
    m4 mtr4 = {{
        { 1.f, 3.f, 4.f, 7.f },
        { 2.f, 6.f, 9.f, 8.f },
        { 5.f, 4.f, 9.f, 1.f },
        { 2.f, 5.f, 2.f, 1.f }
    }};
    m4 rtr4 = m4_transpose(mtr4);
    m4 etr4 = {{
        { 1.f, 2.f, 5.f, 2.f },
        { 3.f, 6.f, 4.f, 5.f },
        { 4.f, 9.f, 9.f, 2.f },
        { 7.f, 8.f, 1.f, 1.f }
    }};
    for_loop_r(4)
        for_loop_c(4)
            assert(rtr4.M[r][c] == etr4.M[r][c]);
    printf("     -> OK!\n");

    printf("\n+ test m*_llt +\n");
    printf(" - test m2 simmetric\n");
    m2 mllt2 = {{
        { 31.f,  9.f },
        {  9.f, 31.f }
    }};
    m2 rllt2 = m2_llt(mllt2);
    m2 tllt2 = m2_mmul(rllt2, m2_transpose(rllt2));
    for_loop_r(2)
        for_loop_c(2)
            assert(tllt2.M[r][c] == mllt2.M[r][c]); // A = L*transpose(L)
    printf("     -> OK!\n");

    printf(" - test m3 simmetric\n");
    m3 mllt3 = {{
        { 31.f,  9.f, 21.f },
        {  9.f, 55.f, 11.f },
        { 21.f, 11.f, 44.f }
    }};
    m3 rllt3 = m3_llt(mllt3);
    m3 tllt3 = m3_mmul(rllt3, m3_transpose(rllt3));
    for_loop_r(3)
        for_loop_c(3)
            assert(tllt3.M[r][c] == mllt3.M[r][c]); // A = L*transpose(L)
    printf("     -> OK!\n");

    printf(" - test m4 simmetric\n");
    m4 mllt4 = {{
        { 18.f, 22.f,  54.f,  42.f },
        { 22.f, 70.f,  86.f,  62.f },
        { 54.f, 86.f, 174.f, 134.f },
        { 42.f, 62.f, 134.f, 106.f  }
    }};
    m4 rllt4 = m4_llt(mllt4);
    m4 tllt4 = m4_mmul(rllt4, m4_transpose(rllt4));
    for_loop_r(4)
        for_loop_c(4)
            assert(fabs(tllt4.M[r][c] - mllt4.M[r][c]) < EPSILON); // A = L*transpose(L)
    printf("     -> OK!\n");

    printf("\n+ test m*s_inv +\n");

    printf(" - test m2\n");
    m2 tbis2 = {{
        { 5.f, 3.f },
        { 3.f, 5.f }
    }};
    m2 invs2 = m2s_inv(tbis2);
    m2 exinvs2 = {{
        {  0.3125f, -0.1875f },
        { -0.1875f,  0.3125f }
    }};
    start_matrix_loop(2)
    {
        assert(fabs(invs2.M[r][c] - exinvs2.M[r][c]) < EPSILON);
    }
    end_matrix_loop
    m2 hopefully_eye2 = m2_mmul(tbis2, invs2);
    m2 eye2_i = m2_eye();
    start_matrix_loop(2)
        assert(fabs(eye2_i.M[r][c] - hopefully_eye2.M[r][c]) < EPSILON);
    end_matrix_loop
    printf("     -> OK!\n");

    printf(" - test m3\n");
    m3 tbis3 = {{
        {  6.f,  15.f,  55.f },
        { 15.f,  55.f, 225.f },
        { 55.f, 225.f, 979.f }
    }};
    m3 invs3 = m3s_inv(tbis3);
    m3 exinvs3 = {{
        {  0.8214285f, -0.5892857f,  0.0892857f },
        { -0.5892857f,  0.7267857f, -0.1339285f },
        {  0.0892857f, -0.1339285f,  0.0267857f }
    }};
    start_matrix_loop(3)
        assert(fabs(invs3.M[r][c] - exinvs3.M[r][c]) < EPSILON);
    end_matrix_loop
    m3 hopefully_eye3 = m3_mmul(tbis3, invs3);
    m3 eye3_i = m3_eye();
    start_matrix_loop(3)
        assert(fabs(eye3_i.M[r][c] - hopefully_eye3.M[r][c]) < EPSILON);
    end_matrix_loop
    printf("     -> OK!\n");

    printf(" - test m4\n");
    m4 tbis4 = {{
        { 18.f, 22.f,  54.f,  42.f },
        { 22.f, 70.f,  86.f,  62.f },
        { 54.f, 86.f, 174.f, 134.f },
        { 42.f, 62.f, 134.f, 106.f }
    }};
    m4 invs4 = m4s_inv(tbis4);
    m4 exinvs4 = {{
        {  2.515625f,  0.484375f, -1.296874f,  0.359374f },
        {  0.484375f,  0.140625f, -0.328125f,  0.140625f },
        { -1.296874f, -0.328125f,  1.015625f, -0.578125f },
        {  0.359374f,  0.140625f, -0.578125f,  0.515625f }
    }};
    start_matrix_loop(4)
        assert(fabs(invs4.M[r][c] - exinvs4.M[r][c]) < EPSILON);
    end_matrix_loop
    m4 hopefully_eye4 = m4_mmul(tbis4, invs4);
    m4 eye4_i = m4_eye();
    start_matrix_loop(4)
    {
        printf("R[%d][%d] -> %f == %f\n", r, c, eye4_i.M[r][c], hopefully_eye4.M[r][c]);
        assert(fabs(eye4_i.M[r][c] - hopefully_eye4.M[r][c]) < EPSILON);
    }
    end_matrix_loop
    printf("     -> OK!\n");

#if 0
    printf("\n+ test m*_inv +\n");
    printf(" - test m2\n");
    printf(" - test m3\n");
    printf(" - test m4\n");
    m4 minv4 = {{
        { 10.f, 210.f, 2.f, 9.f },
        { 0.2f, 33.f, 12.f, 99.f },
        { 2.1f, 99.f, 0.8f, 1.2f },
        { 1.1f, 1.2f, 5.9f, 41.f }
    }};
    m4 ernv4 = {{}};
    m4 rsnv4 = m4_inv(minv4);

    start_matrix_loop(4)
        assert(rsnv4.M[r][c] == ernv4.M[r][c]);
    end_matrix_loop
#endif

    printf("\n+ test v*_ssum +\n");

    printf(" - test v2\n");
    v2 as2 = { 5.f, 6.f };
    v2 rs2 = v2_ssum(2.f, as2);
    assert(rs2.x == 7.f);
    assert(rs2.y == 8.f);
    printf("     -> OK!\n");

    printf(" - test v3\n");
    v3 as3 = { 5.f, 6.f, 7.f };
    v3 rs3 = v3_ssum(2.f, as3);
    v3_ssum(2.f, as3);
    assert(rs3.x == 7.f);
    assert(rs3.y == 8.f);
    assert(rs3.z == 9.f);
    printf("     -> OK!\n");

    printf(" - test v3\n");
    v4 as4 = { 5.f, 6.f, 7.f, 8.f };
    v4 rs4 = v4_ssum(2.f, as4);
    assert(rs4.w == 7.f);
    assert(rs4.x == 8.f);
    assert(rs4.y == 9.f);
    assert(rs4.z == 10.f);
    printf("     -> OK!\n");

    printf("\n+ test m*_ssum +\n");

    printf(" - test m2\n");
    m2 mss2 = {{
        { 10.f, 10.f },
        { 10.f, 10.f }
    }};
    m2 ress2 = m2_ssum(2.f, mss2);
    m2 exrs2 = {{
        { 12.f, 12.f },
        { 12.f, 12.f }
    }};
    for (int c = 0; c < 2; ++c)
    {
        for (int r = 0; r < 2; ++r)
        {
            assert(ress2.M[r][c] == exrs2.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf(" - test m3\n");
    m3 mss3 = {{
        { 10.f, 10.f, 20.f },
        { 10.f, 10.f, 20.f },
        { 10.f, 10.f, 20.f }
    }};
    m3 ress3 = m3_ssum(2.f, mss3);
    m3 exrs3 = {{
        { 12.f, 12.f, 22.f },
        { 12.f, 12.f, 22.f },
        { 12.f, 12.f, 22.f }
    }};
    for (int c = 0; c < 3; ++c)
    {
        for (int r = 0; r < 3; ++r)
        {
            assert(ress3.M[r][c] == exrs3.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf(" - test m4\n");
    m4 mss4 = {{
        { 10.f, 10.f, 20.f, 100.f },
        { 10.f, 10.f, 20.f, 100.f  },
        { 10.f, 10.f, 20.f, 100.f },
        { 10.f, 10.f, 20.f, 100.f }
    }};
    m4 ress4 = m4_ssum(2.f, mss4);
    m4 exrs4 = {{
        { 12.f, 12.f, 22.f, 102.f },
        { 12.f, 12.f, 22.f, 102.f  },
        { 12.f, 12.f, 22.f, 102.f },
        { 12.f, 12.f, 22.f, 102.f }
    }};
    for (int c = 0; c < 4; ++c)
    {
        for (int r = 0; r < 4; ++r)
        {
            assert(ress4.M[r][c] == exrs4.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf("\n+ test v*_msum +\n");

    printf(" - test v2\n");
    v2 amsv2 = { 1.f, 5.f };
    v2 bmsv2 = { 3.f, 2.f };
    v2 rmsv2 = v2_msum(amsv2, bmsv2);
    assert(rmsv2.x == amsv2.x + bmsv2.x);
    assert(rmsv2.y == amsv2.y + bmsv2.y);
    printf("     -> OK!\n");

    printf(" - test v3\n");
    v3 amsv3 = { 1.f, 5.f, 7.f };
    v3 bmsv3 = { 3.f, 2.f, 9.f };
    v3 rmsv3 = v3_msum(amsv3, bmsv3);
    assert(rmsv3.x == amsv3.x + bmsv3.x);
    assert(rmsv3.y == amsv3.y + bmsv3.y);
    assert(rmsv3.z == amsv3.z + bmsv3.z);
    printf("     -> OK!\n");

    printf(" - test v4\n");
    v4 amsv4 = { 1.f, 5.f, 7.f, 1.f };
    v4 bmsv4 = { 3.f, 2.f, 9.f, 5.f };
    v4 rmsv4 = v4_msum(amsv4, bmsv4);
    assert(rmsv4.w == amsv4.w + bmsv4.w);
    assert(rmsv4.x == amsv4.x + bmsv4.x);
    assert(rmsv4.y == amsv4.y + bmsv4.y);
    assert(rmsv4.z == amsv4.z + bmsv4.z);
    printf("     -> OK!\n");

    printf("\n+ test m*_msum +\n");
    printf(" - test m2\n");
    m2 amsm2 = {{
        { 1.f, 2.f },
        { 1.f, 2.f }
    }};
    m2 bmsm2 = {{
        { 5.f, 6.f },
        { 5.f, 6.f }
    }};
    m2 rmsm2 = m2_msum(amsm2, bmsm2);
    for (int c = 0; c < 2; ++c)
    {
        for (int r = 0; r < 2; ++r)
        {
            assert(rmsm2.M[r][c] == amsm2.M[r][c] + bmsm2.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf(" - test m3\n");
    m3 amsm3 = {{
        { 1.f, 2.f, 3.f },
        { 1.f, 2.f, 3.f },
        { 1.f, 2.f, 3.f }
    }};
    m3 bmsm3 = {{
        { 5.f, 6.f, 6.f },
        { 5.f, 6.f, 6.f },
        { 5.f, 6.f, 6.f }
    }};
    m3 rmsm3 = m3_msum(amsm3, bmsm3);
    for (int c = 0; c < 3; ++c)
    {
        for (int r = 0; r < 3; ++r)
        {
            assert(rmsm3.M[r][c] == amsm3.M[r][c] + bmsm3.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf(" - test m4\n");
    m4 amsm4 = {{
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f }
    }};
    m4 bmsm4 = {{
        { 5.f, 6.f, 6.f, 7.f },
        { 5.f, 6.f, 6.f, 7.f },
        { 5.f, 6.f, 6.f, 7.f },
        { 5.f, 6.f, 6.f, 7.f }
    }};
    m4 rmsm4 = m4_msum(amsm4, bmsm4);
    for (int c = 0; c < 4; ++c)
    {
        for (int r = 0; r < 4; ++r)
        {
            assert(rmsm4.M[r][c] == amsm4.M[r][c] + bmsm4.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf("\n+ test v*_smul +\n");

    printf(" - test v2\n");
    v2 vsmul2 = { 1.f, 5.f };
    v2 rsmul2 = v2_smul(2.f, vsmul2);
    assert(rsmul2.x == 2.f * vsmul2.x);
    assert(rsmul2.y == 2.f * vsmul2.y);
    printf("     -> OK!\n");

    printf(" - test v3\n");
    v3 vsmul3 = { 1.f, 5.f, 7.f };
    v3 rsmul3 = v3_smul(2.f, vsmul3);
    assert(rsmul3.x == 2.f * vsmul3.x);
    assert(rsmul3.y == 2.f * vsmul3.y);
    assert(rsmul3.z == 2.f * vsmul3.z);
    printf("     -> OK!\n");

    printf(" - test v4\n");
    v4 vsmul4 = { 1.f, 5.f, 7.f, 1.f };
    v4 rsmul4 = v4_smul(2.f, vsmul4);
    assert(rsmul4.w == 2.f * vsmul4.w);
    assert(rsmul4.x == 2.f * vsmul4.x);
    assert(rsmul4.y == 2.f * vsmul4.y);
    assert(rsmul4.z == 2.f * vsmul4.z);
    printf("     -> OK!\n");

    printf("\n+ test m*_smul +\n");
    
    printf(" - test m2\n");
    m2 msmu2 = {{
        { 1.f, 2.f },
        { 1.f, 2.f }
    }};
    m2 rsmu2 = m2_smul(2.f, msmu2);
    for (int c = 0; c < 2; ++c)
    {
        for (int r = 0; r < 2; ++r)
        {
            assert(rsmu2.M[r][c] == 2.f * msmu2.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf(" - test v3\n");
    m3 msmu3 = {{
        { 1.f, 2.f, 3.f },
        { 1.f, 2.f, 3.f },
        { 1.f, 2.f, 3.f }
    }};
    m3 rsmu3 = m3_smul(2.f, msmu3);
    for (int c = 0; c < 3; ++c)
    {
        for (int r = 0; r < 3; ++r)
        {
            assert(rsmu3.M[r][c] == 2.f * msmu3.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf(" - test v4\n");
    m4 msmu4 = {{
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f }
    }};
    m4 rsmu4 = m4_smul(2.f, msmu4);
    for (int c = 0; c < 4; ++c)
    {
        for (int r = 0; r < 4; ++r)
        {
            assert(rsmu4.M[r][c] == 2.f * msmu4.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf("\n+ test *2_mmul +\n");

    printf(" - test m2\n");
    m2 ammul2 = {{
        { 1.f, 4.f },
        { 3.f, 2.f }
    }};
    m2 bmmul2 = {{
        { 5.f,  9.f },
        { 6.f, 10.f }
    }};
    m2 ermmul2 = {{
        { 29.f, 49.f },
        { 27.f, 47.f }
    }};
    m2 remmul2 = m2_mmul(ammul2, bmmul2); // 2x2, 2x2 -> 2x2
    for (int c = 0; c < 2; ++c)
    {
        for (int r = 0; r < 2; ++r)
        {
            assert(remmul2.M[r][c] == ermmul2.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf(" - test v2\n");
    v2 avmul2 = { 1.f, 2.f };
    v2 bvmul2 = { 5.f, 9.f };
    m2 ervmul2 = {{
        {  5.f,  9.f },
        { 10.f, 18.f }
    }};
    m2 revmul2 = v2_mmul(avmul2, bvmul2); // 2x1, 1x2 -> 2x2
    for (int c = 0; c < 2; ++c)
    {
        for (int r = 0; r < 2; ++r)
        {
            assert(revmul2.M[r][c] == ervmul2.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf("\n+ test *2*2_mmul +\n");
    m2 mmul2 = {{
        { 1.f, 2.f },
        { 1.f, 2.f }
    }};
    v2 vmul2 = { 12.f, 4.f };

    printf(" - test v2m2\n");
    v2 ervmmul2 = { 16.f, 32.f };
    v2 rsvmmul2 = v2m2_mmul(vmul2, mmul2); // 1x2, 2x2 -> 1x2
    assert(ervmmul2.x == rsvmmul2.x);
    assert(ervmmul2.y == rsvmmul2.y);
    printf("     -> OK!\n");

    printf(" - test m2v2\n");
    v2 ermvmul2 = { 20.f, 20.f };
    v2 rsmvmul2 = m2v2_mmul(mmul2, vmul2); // 2x2, 2x1 -> 2x1
    assert(ermvmul2.x == rsmvmul2.x);
    assert(ermvmul2.y == rsmvmul2.y);
    printf("     -> OK!\n");

    printf("\n+ test *3_mmul +\n");

    printf(" - test m3\n");
    m3 ammul3 = {{
        { 1.f, 2.f, 3.f },
        { 4.f, 6.f, 8.f },
        { 5.f, 7.f, 9.f }   
    }};
    m3 bmmul3 = {{
        { 5.f,  9.f, 13.f },
        { 6.f, 10.f, 14.f },
        { 7.f, 11.f, 15.f }
    }};
    m3 ermmul3 = {{
        {  38.f,  62.f,  86.f },
        { 112.f, 184.f, 256.f },
        { 130.f, 214.f, 298.f }
    }};
    m3 remmul3 = m3_mmul(ammul3, bmmul3); // 3x3, 3x3 -> 3x3
    for (int c = 0; c < 3; ++c)
    {
        for (int r = 0; r < 3; ++r)
        {
            assert(remmul3.M[r][c] == ermmul3.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf(" - test v3\n");
    v3 avmul3 = { 1.f, 2.f,  3.f };
    v3 bvmul3 = { 5.f, 9.f, 13.f };
    m3 ervmul3 = {{
        {  5.f,  9.f, 13.f },
        { 10.f, 18.f, 26.f },
        { 15.f, 27.f, 39.f }
    }};
    m3 revmul3 = v3_mmul(avmul3, bvmul3); // 3x1, 1x3 -> 3x3
    for (int c = 0; c < 3; ++c)
    {
        for (int r = 0; r < 3; ++r)
        {
            assert(revmul3.M[r][c] == ervmul3.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf("\n+ test *3*3_mmul +\n");
    printf(" - test v3m3\n");
    m3 mmul3 = {{
        { 1.f, 2.f, 3.f },
        { 1.f, 2.f, 3.f },
        { 1.f, 2.f, 3.f }
    }};
    v3 vmul3 = { 12.f, 4.f, 3.f };

    v3 ervmmul3 = { 19.f, 38.f, 57.f };
    v3 rsvmmul3 = v3m3_mmul(vmul3, mmul3); // 1x3, 3x3 -> 1x3
    assert(ervmmul3.x == rsvmmul3.x);
    assert(ervmmul3.y == rsvmmul3.y);
    assert(ervmmul3.z == rsvmmul3.z);
    printf("     -> OK!\n");

    printf(" - test m3v3\n");
    v3 ermvmul3 = { 29.f, 29.f, 29.f };
    v3 rsmvmul3 = m3v3_mmul(mmul3, vmul3); // 4x4, 4x1 -> 4x1
    assert(ermvmul3.x == rsmvmul3.x);
    assert(ermvmul3.y == rsmvmul3.y);
    assert(ermvmul3.z == rsmvmul3.z);
    printf("     -> OK!\n");

    printf("\n+ test *4_mmul +\n");

    printf(" - test m4\n");
    m4 ammul4 = {{
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f }
    }};
    m4 bmmul4 = {{
        { 5.f,  9.f, 13.f, 17.f },
        { 6.f, 10.f, 14.f, 18.f },
        { 7.f, 11.f, 15.f, 19.f },
        { 8.f, 12.f, 16.f, 20.f }
    }};
    m4 ermmul4 = {{
        { 70.f, 110.f, 150.f, 190.f },
        { 70.f, 110.f, 150.f, 190.f },
        { 70.f, 110.f, 150.f, 190.f },
        { 70.f, 110.f, 150.f, 190.f }
    }};
    m4 remmul4 = m4_mmul(ammul4, bmmul4); // 4x4, 4x4 -> 4x4
    for (int c = 0; c < 4; ++c)
    {
        for (int r = 0; r < 4; ++r)
        {
            assert(remmul4.M[r][c] == ermmul4.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf(" - test v4\n");
    v4 avmul4 = { 1.f, 2.f, 3.f, 4.f };
    v4 bvmul4 = { 5.f,  9.f, 13.f, 17.f };
    m4 ervmul4 = {{
        {  5.f,  9.f, 13.f, 17.f },
        { 10.f, 18.f, 26.f, 34.f },
        { 15.f, 27.f, 39.f, 51.f },
        { 20.f, 36.f, 52.f, 68.f }
    }};
    m4 revmul4 = v4_mmul(avmul4, bvmul4); // 4x1, 1x4 -> 4x4
    for (int c = 0; c < 4; ++c)
    {
        for (int r = 0; r < 4; ++r)
        {
            assert(revmul4.M[r][c] == ervmul4.M[r][c]);
        }
    }
    printf("     -> OK!\n");

    printf("\n+ test *4*4_mmul +\n");
    printf(" - test v4m4\n");
    m4 mmul4 = {{
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f },
        { 1.f, 2.f, 3.f, 4.f }
    }};
    v4 vmul4 = { 12.f, 4.f, 3.f, 0.6f };
    v4 ervmmul4 = { 19.6f, 39.2f, 58.8f, 78.4f };
    v4 rsvmmul4 = v4m4_mmul(vmul4, mmul4); // 1x4, 4x4 -> 1x4
    assert(ervmmul4.w == rsvmmul4.w);
    assert(ervmmul4.x == rsvmmul4.x);
    assert(ervmmul4.y == rsvmmul4.y);
    assert(ervmmul4.z == rsvmmul4.z);
    printf("     -> OK!\n");

    printf(" - test m4v4\n");
    v4 ermvmul4 = { 31.4f, 31.4f, 31.4f, 31.4f };
    v4 rsmvmul4 = m4v4_mmul(mmul4, vmul4); // 4x4, 4x1 -> 4x1
    assert(ermvmul4.w == rsmvmul4.w);
    assert(ermvmul4.x == rsmvmul4.x);
    assert(ermvmul4.y == rsmvmul4.y);
    assert(ermvmul4.z == rsmvmul4.z);
    printf("     -> OK!\n");
}

void test_afrewrls()
{
    afewrls ls;
    f32 alpha = 0.95;
    f32 lambda_minus = 0.7f;
    f32 lambda_plus = 1.f;
    afewrls_init(&ls, alpha, lambda_minus, lambda_plus);
    print_afewrls(&ls);

    v4 x1 = { 10.f, 12.f, 14.f, 3.f };

    afewrls_update(&ls, x1, 12.f);
    print_afewrls(&ls);

    v4 x2 = { 10.1f, 11.8f, 9.f, 3.5f };
    afewrls_update(&ls, x2, 12.f);
    print_afewrls(&ls);
    afewrls_update(&ls, x2, 11.f);
    print_afewrls(&ls);

    v4 x_ = {1.f, 1.f, 1.f, 1.f};
    afewrls_predict(&ls, x_);
#define I 100
    for (int i = 0; i < I; ++i)
    {
        afewrls_update(&ls, x_, 11.f);
    }
    f32 res = afewrls_predict(&ls, x_);
    printf("res: %.5f\n", res);

}

void test_fls()
{
    printf("\n\ntesting fls:\n");
    fls ls;
    f32 factor = 0.95f;
    fls_init(&ls, factor);

    v4 x1 = { 10.f, 12.f, 14.f, 3.f };
    for (int i = 0; i < 50; ++i)
    {
        printf("i: %d\n", i);
        fls_update(&ls, x1, 12.0f);
    }

    print4(ls.beta);
}

int main()
{
    test_la();
    // test_afrewrls();
    test_fls();
}