
#ifndef FASTTREE_SSE128OPERATIONS_TCC
#define FASTTREE_SSE128OPERATIONS_TCC

#include "SSE128Operations.h"

template<>
template<>
inline float fasttree::SSE128Operations<float>::mm_sum(register __m128 sum) {
    /* stupider but faster */
    alignas(ALIGNMENT) float f[4];
    _mm_store_ps(f, sum);
    return (f[0] + f[1] + f[2] + f[3]);
}

template<>
template<>
inline double fasttree::SSE128Operations<double>::mm_sum(register __m128d sum) {
    /* stupider but faster */
    alignas(ALIGNMENT) double f[2];;
    _mm_store_pd(f, sum);
    return (f[0] + f[1]);
}

template<>
void fasttree::SSE128Operations<float>::
vector_multiply(float f1[], float f2[], int64_t n, float fOut[]) {
    for (int64_t i = 0; i < n; i += 4) {
        __m128 a, b, c;
        a = _mm_load_ps(f1 + i);
        b = _mm_load_ps(f2 + i);
        c = _mm_mul_ps(a, b);
        _mm_store_ps(fOut + i, c);
    }
}

template<>
void fasttree::SSE128Operations<double>::
vector_multiply(double f1[], double f2[], int64_t n, double fOut[]) {
    for (int64_t i = 0; i < n; i += 2) {
        __m128d a, b, c;
        a = _mm_load_pd(f1 + i);
        b = _mm_load_pd(f2 + i);
        c = _mm_mul_pd(a, b);
        _mm_store_pd(fOut + i, c);
    }
}

template<>
float fasttree::SSE128Operations<float>
::vector_multiply_sum(float f1[], float f2[], int64_t n) {
    if (n == 4) {
        return f1[0] * f2[0] + f1[1] * f2[1] + f1[2] * f2[2] + f1[3] * f2[3];
    }
    __m128 sum = _mm_setzero_ps();
    for (int64_t i = 0; i < n; i += 4) {
        __m128 a, b, c;
        a = _mm_load_ps(f1 + i);
        b = _mm_load_ps(f2 + i);
        c = _mm_mul_ps(a, b);
        sum = _mm_add_ps(c, sum);
    }
    return mm_sum(sum);
}

template<>
double fasttree::SSE128Operations<double>
::vector_multiply_sum(double f1[], double f2[], int64_t n) {
    if (n == 4) {
        return f1[0] * f2[0] + f1[1] * f2[1] + f1[2] * f2[2] + f1[3] * f2[3];
    }
    __m128d sum = _mm_setzero_pd();
    for (int64_t i = 0; i < n; i += 2) {
        __m128d a, b, c;
        a = _mm_load_pd(f1 + i);
        b = _mm_load_pd(f2 + i);
        c = _mm_mul_pd(a, b);
        sum = _mm_add_pd(c, sum);
    }
    return mm_sum(sum);
}

template<>
float fasttree::SSE128Operations<float>::
vector_multiply3_sum(float f1[], float f2[], float f3[], int64_t n) {
    __m128 sum = _mm_setzero_ps();
    for (int64_t i = 0; i < n; i += 4) {
        __m128 a1, a2, a3;
        a1 = _mm_load_ps(f1 + i);
        a2 = _mm_load_ps(f2 + i);
        a3 = _mm_load_ps(f3 + i);
        sum = _mm_add_ps(_mm_mul_ps(_mm_mul_ps(a1, a2), a3), sum);
    }
    return mm_sum(sum);
}

template<>
double fasttree::SSE128Operations<double>::
vector_multiply3_sum(double f1[], double f2[], double f3[], int64_t n) {
    __m128d sum = _mm_setzero_pd();
    for (int64_t i = 0; i < n; i += 2) {
        __m128d a1, a2, a3;
        a1 = _mm_load_pd(f1 + i);
        a2 = _mm_load_pd(f2 + i);
        a3 = _mm_load_pd(f3 + i);
        sum = _mm_add_pd(_mm_mul_pd(_mm_mul_pd(a1, a2), a3), sum);
    }
    return mm_sum(sum);
}

template<>
float fasttree::SSE128Operations<float>::
vector_dot_product_rot(float f1[], float f2[], float fBy[], int64_t n) {
    __m128 sum1 = _mm_setzero_ps();
    __m128 sum2 = _mm_setzero_ps();
    for (int64_t i = 0; i < n; i += 4) {
        __m128 a1, a2, aBy;
        a1 = _mm_load_ps(f1 + i);
        a2 = _mm_load_ps(f2 + i);
        aBy = _mm_load_ps(fBy + i);
        sum1 = _mm_add_ps(_mm_mul_ps(a1, aBy), sum1);
        sum2 = _mm_add_ps(_mm_mul_ps(a2, aBy), sum2);
    }
    return (mm_sum(sum1) * mm_sum(sum2));
}

template<>
double fasttree::SSE128Operations<double>::
vector_dot_product_rot(double f1[], double f2[], double fBy[], int64_t n) {
    __m128d sum1 = _mm_setzero_pd();
    __m128d sum2 = _mm_setzero_pd();
    for (int64_t i = 0; i < n; i += 2) {
        __m128d a1, a2, aBy;
        a1 = _mm_load_pd(f1 + i);
        a2 = _mm_load_pd(f2 + i);
        aBy = _mm_load_pd(fBy + i);
        sum1 = _mm_add_pd(_mm_mul_pd(a1, aBy), sum1);
        sum2 = _mm_add_pd(_mm_mul_pd(a2, aBy), sum2);
    }
    return mm_sum(sum1) * mm_sum(sum2);
}

template<>
float fasttree::SSE128Operations<float>::
vector_sum(float f1[], int64_t n) {
    if (n == 4) {
        return f1[0] + f1[1] + f1[2] + f1[3];
    }
    __m128 sum = _mm_setzero_ps();
    for (int64_t i = 0; i < n; i += 4) {
        __m128 a;
        a = _mm_load_ps(f1 + i);
        sum = _mm_add_ps(a, sum);
    }
    return (mm_sum(sum));
}

template<>
double fasttree::SSE128Operations<double>::
vector_sum(double f1[], int64_t n) {
    if (n == 4) {
        return f1[0] + f1[1] + f1[2] + f1[3];
    }
    __m128d sum = _mm_setzero_pd();
    for (int64_t i = 0; i < n; i += 2) {
        __m128d a;
        a = _mm_load_pd(f1 + i);
        sum = _mm_add_pd(a, sum);
    }
    return mm_sum(sum);
}

template<>
void fasttree::SSE128Operations<float>::
vector_multiply_by(float f[], float fBy, int64_t n, float fOut[]) {
    __m128 c = _mm_set1_ps(fBy);
    for (int64_t i = 0; i < n; i += 4) {
        __m128 a, b;
        a = _mm_load_ps(f + i);
        b = _mm_mul_ps(a, c);
        _mm_store_ps(fOut + i, b);
    }
}

template<>
void fasttree::SSE128Operations<double>::
vector_multiply_by(double f[], double fBy, int64_t n, double fOut[]) {
    __m128d c = _mm_set1_pd(fBy);
    for (int64_t i = 0; i < n; i += 2) {
        __m128d a, b;
        a = _mm_load_pd(f + i);
        b = _mm_mul_pd(a, c);
        _mm_store_pd(fOut + i, b);
    }
}

template<>
void fasttree::SSE128Operations<float>::
vector_add_mult(float fTot[], float fAdd[], float weight, int64_t n) {
    __m128 w = _mm_set1_ps(weight);
    for (int64_t i = 0; i < n; i += 4) {
        __m128 tot, add;
        tot = _mm_load_ps(fTot + i);
        add = _mm_load_ps(fAdd + i);
        _mm_store_ps(fTot + i, _mm_add_ps(tot, _mm_mul_ps(add, w)));
    }
}

template<>
void fasttree::SSE128Operations<double>::
vector_add_mult(double fTot[], double fAdd[], double weight, int64_t n) {
    __m128d w = _mm_set1_pd(weight);
    for (int64_t i = 0; i < n; i += 2) {
        __m128d tot, add;
        tot = _mm_load_pd(fTot + i);
        add = _mm_load_pd(fAdd + i);
        _mm_store_pd(fTot + i, _mm_add_pd(tot, _mm_mul_pd(add, w)));
    }
}

template<>
void fasttree::SSE128Operations<float>::
matrixt_by_vector4(float mat[4][4], float vec[4], float out[4]) {
    __m128 o = _mm_setzero_ps();
    /* result is a sum of vectors: sum(k) v[k] * mat[k][] */
    for (int64_t j = 0; j < 4; j++) {
        __m128 m = _mm_load_ps(&mat[j][0]);
        __m128 vj = _mm_load1_ps(&vec[j]);    /* is it faster to shuffle v? */
        o = _mm_add_ps(o, _mm_mul_ps(vj, m));
    }
    _mm_store_ps(out, o);
}

template<>
void fasttree::SSE128Operations<double>::
matrixt_by_vector4(double mat[4][4], double vec[4], double out[4]) {
    __m128d o1 = _mm_setzero_pd();
    __m128d o2 = _mm_setzero_pd();
    /* result is a sum of vectors: sum(k) v[k] * mat[k][] */
    for (int64_t j = 0; j < 4; j++) {
        __m128d vj = _mm_load1_pd(&vec[j]); /* is it faster to shuffle v? */
        __m128d m1 = _mm_load_pd(&mat[j][0]);
        __m128d m2 = _mm_load_pd(&mat[j][2]);
        o1 = _mm_add_pd(o1, _mm_mul_pd(vj, m1));
        o2 = _mm_add_pd(o2, _mm_mul_pd(vj, m2));
    }
    _mm_store_pd(out, o1);
    _mm_store_pd(&out[2], o2);
}

template<>
void fasttree::SSE128Operations<float>::fastexp(float fTot[], int64_t n, int lvl) {
    if (lvl == 0) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = (float) std::exp((double) fTot[k]);
        }
    } else if (lvl == 1) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = std::exp(fTot[k]);
        }
    } else if (lvl == 2) {
        alignas(ALIGNMENT) double arr[2];
        for (int64_t k = 0; k < n; k += 2) {
            __m128d nums = _mm_set_pd(fTot[k + 1], fTot[k]);
            __m128d res = fastexpImpl(nums);
            _mm_store_pd(arr, res);
            fTot[k] = (float) arr[0];
            fTot[k + 1] = (float) arr[1];
        }
    } else {
        for (int64_t k = 0; k < n; k += 4) {
            __m128 nums = _mm_load_ps(&fTot[k]);
            __m128 res = fastexpImpl(nums);
            _mm_store_ps(&fTot[k], res);
        }
    }
}

template<>
void fasttree::SSE128Operations<double>::fastexp(double fTot[], int64_t n, int lvl) {
    if (lvl == 0) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = (numeric_t) std::exp((double) fTot[k]);
        }
    } else if (lvl == 1) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = (numeric_t) std::exp((float) fTot[k]);
        }
    } else if (lvl == 2) {
        for (int64_t k = 0; k < n; k += 2) {
            __m128d nums = _mm_load_pd(&fTot[k]);
            __m128d res = fastexpImpl(nums);
            _mm_store_pd(&fTot[k], res);
        }
    } else {
        alignas(ALIGNMENT) float arr[4];
        for (int64_t k = 0; k < n; k += 4) {
            __m128 nums = _mm_set_ps((float) fTot[k + 3], (float) fTot[k + 2], (float) fTot[k + 1], (float) fTot[k]);
            __m128 res = fastexpImpl(nums);
            _mm_store_ps(arr, res);
            fTot[k] = (float) arr[0];
            fTot[k + 1] = (float) arr[1];
            fTot[k + 2] = (float) arr[2];
            fTot[k + 3] = (float) arr[3];
        }
    }
}

template<typename Precision>
__m128 fasttree::SSE128Operations<Precision>::fastexpImpl(__m128 x) {
    __m128i m, u;
    __m128 xx, px, qx;

    px = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(1.4426950408889634073599f), x), _mm_set1_ps(0.5));
    #if (defined __SSE4_1__) || (defined __AVX__)
    px = _mm_floor_ps(px);
    #else
    __m128 j = _mm_set1_ps(1.0f);
    __m128i i = _mm_cvttps_epi32(px);
    __m128 fi = _mm_cvtepi32_ps(i);
    __m128 igx = _mm_cmpgt_ps(fi, px);
    j = _mm_and_ps(igx, j);
    px = _mm_sub_ps(fi, j);
    #endif

    /* m = (__m128i) px */
    m = _mm_cvttps_epi32(px);

    x = _mm_sub_ps(x, _mm_mul_ps(px, _mm_set1_ps(6.93145751953125E-1f)));
    x = _mm_sub_ps(x, _mm_mul_ps(px, _mm_set1_ps(1.42860682030941723212E-6f)));
    xx = _mm_mul_ps(x, x);

    /* px = x * P(x**2). */
    px = _mm_set1_ps(1.26177193074810590878E-4f);
    px = _mm_mul_ps(px, xx);
    px = _mm_add_ps(px, _mm_set1_ps(3.02994407707441961300E-2f));
    px = _mm_mul_ps(px, xx);
    px = _mm_add_ps(px, _mm_set1_ps(9.99999999999999999910E-1f));
    px = _mm_mul_ps(px, x);

    /* Evaluate Q(x**2). */
    qx = _mm_set1_ps(3.00198505138664455042E-6f);
    qx = _mm_mul_ps(qx, xx);
    qx = _mm_add_ps(qx, _mm_set1_ps(2.52448340349684104192E-3f));
    qx = _mm_mul_ps(qx, xx);
    qx = _mm_add_ps(qx, _mm_set1_ps(2.27265548208155028766E-1f));
    qx = _mm_mul_ps(qx, xx);
    qx = _mm_add_ps(qx, _mm_set1_ps(2.00000000000000000009E0f));

    /* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) ) */
    x = _mm_div_ps(px, _mm_sub_ps(qx, px));
    x = _mm_add_ps(_mm_set1_ps(1.0), _mm_mul_ps(_mm_set1_ps(2.0), x));

    /* Build 2^n in double. */
    m = _mm_add_epi32(m, _mm_set1_epi32(127));
    u = _mm_slli_epi32(m, 23);

    return _mm_mul_ps(x, _mm_castsi128_ps(u));
}

template<typename Precision>
__m128d fasttree::SSE128Operations<Precision>::fastexpImpl(__m128d x) {
    __m128i m, u;
    __m128d xx, px, qx;

    px = _mm_add_pd(_mm_mul_pd(_mm_set1_pd(1.4426950408889634073599), x), _mm_set1_pd(0.5));
    #if (defined __SSE4_1__) || (defined __AVX__)
    px = _mm_floor_pd(px);
    #else
    __m128d j = _mm_set1_pd(1.0f);

    __m128i i = _mm_sub_epi64(
            _mm_castpd_si128(_mm_add_pd(px, _mm_set1_pd(0x0018000000000000))),
            _mm_castpd_si128(_mm_set1_pd(0x0018000000000000))
    );

    __m128d fi = _mm_sub_pd(
            _mm_castsi128_pd(_mm_add_epi64(i, _mm_castpd_si128(_mm_set1_pd(0x0018000000000000)))),
            _mm_set1_pd(0x0018000000000000)
    );

    __m128d igx = _mm_cmpgt_pd(fi, px);
    j = _mm_and_pd(igx, j);
    px = _mm_sub_pd(fi, j);
    #endif

    /* m = (__m128i) px */
    m = _mm_sub_epi64(
            _mm_castpd_si128(_mm_add_pd(px, _mm_set1_pd(0x0018000000000000))),
            _mm_castpd_si128(_mm_set1_pd(0x0018000000000000))
    );

    x = _mm_sub_pd(x, _mm_mul_pd(px, _mm_set1_pd(6.93145751953125E-1)));
    x = _mm_sub_pd(x, _mm_mul_pd(px, _mm_set1_pd(1.42860682030941723212E-6)));
    xx = _mm_mul_pd(x, x);

    /* px = x * P(x**2). */
    px = _mm_set1_pd(1.26177193074810590878E-4);
    px = _mm_mul_pd(px, xx);
    px = _mm_add_pd(px, _mm_set1_pd(3.02994407707441961300E-2));
    px = _mm_mul_pd(px, xx);
    px = _mm_add_pd(px, _mm_set1_pd(9.99999999999999999910E-1));
    px = _mm_mul_pd(px, x);

    /* Evaluate Q(x**2). */
    qx = _mm_set1_pd(3.00198505138664455042E-6);
    qx = _mm_mul_pd(qx, xx);
    qx = _mm_add_pd(qx, _mm_set1_pd(2.52448340349684104192E-3));
    qx = _mm_mul_pd(qx, xx);
    qx = _mm_add_pd(qx, _mm_set1_pd(2.27265548208155028766E-1));
    qx = _mm_mul_pd(qx, xx);
    qx = _mm_add_pd(qx, _mm_set1_pd(2.00000000000000000009E0));

    /* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) ) */
    x = _mm_div_pd(px, _mm_sub_pd(qx, px));
    x = _mm_add_pd(_mm_set1_pd(1.0), _mm_mul_pd(_mm_set1_pd(2.0), x));

    /* Build 2^n in double. */
    m = _mm_add_epi64(m, _mm_set1_epi64x(1023));
    u = _mm_slli_epi64(m, 52);

    return _mm_mul_pd(x, _mm_castsi128_pd(u));
}


#endif

