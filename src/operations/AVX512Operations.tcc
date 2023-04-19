
#include "AVX512Operations.h"
#include <cmath>

#if (defined _WIN32 || defined WIN32 || defined WIN64 || defined _WIN64) //Fix windows union operator[]
#define getf(array, i) ((float*)&array)[i]
#define getd(array, i) ((double*)&array)[i]
#else
#define getf(array, i) array[i]
#define getd(array, i) array[i]
#endif

template<>
inline float veryfasttree::AVX512Operations<float>::mm_sum(__m128 sum) {
    sum = _mm_hadd_ps(sum, sum);
    return getf(sum, 0) + getf(sum, 1);
}

template<>
inline double veryfasttree::AVX512Operations<double>::mm_sum(__m256d sum) {
    sum = _mm256_hadd_pd(sum, sum);
    return getd(sum, 0) + getd(sum, 2);
}

template<>
inline void
veryfasttree::AVX512Operations<float>::vector_multiply(float f1[], float f2[], int64_t n, float fOut[]) {
    int64_t m = n - (n % 16);
    for (int64_t i = 0; i < m; i += 16) {
        __m512 a, b, c;
        a = _mm512_load_ps(f1 + i);
        b = _mm512_load_ps(f2 + i);
        c = _mm512_mul_ps(a, b);
        _mm512_store_ps(fOut + i, c);
    }
    __m128 a, b, c;
    a = _mm_load_ps(f1 + m);
    b = _mm_load_ps(f2 + m);
    c = _mm_mul_ps(a, b);
    _mm_store_ps(fOut + m, c);
}

template<>
inline void
veryfasttree::AVX512Operations<double>::vector_multiply(double f1[], double f2[], int64_t n, double fOut[]) {
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m512d a, b, c;
        a = _mm512_load_pd(f1 + i);
        b = _mm512_load_pd(f2 + i);
        c = _mm512_mul_pd(a, b);
        _mm512_store_pd(fOut + i, c);
    }
    __m256d a, b, c;
    a = _mm256_load_pd(f1 + m);
    b = _mm256_load_pd(f2 + m);
    c = _mm256_mul_pd(a, b);
    _mm256_store_pd(fOut + m, c);
}

template<>
inline float veryfasttree::AVX512Operations<float>::vector_multiply_sum(float f1[], float f2[], int64_t n) {
    __m512 sum = _mm512_setzero_ps();
    int64_t m = n - (n % 16);
    for (int64_t i = 0; i < m; i += 16) {
        __m512 a, b, c;
        a = _mm512_load_ps(f1 + i);
        b = _mm512_load_ps(f2 + i);
        c = _mm512_mul_ps(a, b);
        sum = _mm512_add_ps(c, sum);
    }
    __m128 a, b, c;
    a = _mm_load_ps(f1 + m);
    b = _mm_load_ps(f2 + m);
    c = _mm_mul_ps(a, b);
    if (n == 4) {
        return mm_sum(c);
    }
    return mm_sum(c) + _mm512_reduce_add_ps(sum);
}

template<>
inline double veryfasttree::AVX512Operations<double>::vector_multiply_sum(double f1[], double f2[], int64_t n) {
    __m512d sum = _mm512_setzero_pd();
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m512d a, b, c;
        a = _mm512_load_pd(f1 + i);
        b = _mm512_load_pd(f2 + i);
        c = _mm512_mul_pd(a, b);
        sum = _mm512_add_pd(c, sum);
    }
    __m256d a, b, c;
    a = _mm256_load_pd(f1 + m);
    b = _mm256_load_pd(f2 + m);
    c = _mm256_mul_pd(a, b);
    if (n == 4) {
        return mm_sum(c);
    }
    return mm_sum(c) + _mm512_reduce_add_pd(sum);
}

template<>
inline float
veryfasttree::AVX512Operations<float>::vector_multiply3_sum(float f1[], float f2[], float f3[], int64_t n) {
    __m512 sum = _mm512_setzero_ps();
    int64_t m = n - (n % 16);
    for (int64_t i = 0; i < m; i += 16) {
        __m512 a1, a2, a3;
        a1 = _mm512_load_ps(f1 + i);
        a2 = _mm512_load_ps(f2 + i);
        a3 = _mm512_load_ps(f3 + i);
        sum = _mm512_add_ps(_mm512_mul_ps(_mm512_mul_ps(a1, a2), a3), sum);
    }
    __m128 a1, a2, a3, r;
    a1 = _mm_load_ps(f1 + m);
    a2 = _mm_load_ps(f2 + m);
    a3 = _mm_load_ps(f3 + m);
    r = _mm_mul_ps(_mm_mul_ps(a1, a2), a3);
    if (n == 4) {
        return mm_sum(r);
    }
    return mm_sum(r) + _mm512_reduce_add_ps(sum);
}

template<>
inline double
veryfasttree::AVX512Operations<double>::vector_multiply3_sum(double f1[], double f2[], double f3[], int64_t n) {
    __m512d sum = _mm512_setzero_pd();
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m512d a1, a2, a3;
        a1 = _mm512_load_pd(f1 + i);
        a2 = _mm512_load_pd(f2 + i);
        a3 = _mm512_load_pd(f3 + i);
        sum = _mm512_add_pd(_mm512_mul_pd(_mm512_mul_pd(a1, a2), a3), sum);
    }
    __m256d a1, a2, a3, r;
    a1 = _mm256_load_pd(f1 + m);
    a2 = _mm256_load_pd(f2 + m);
    a3 = _mm256_load_pd(f3 + m);
    r = _mm256_mul_pd(_mm256_mul_pd(a1, a2), a3);
    if (n == 4) {
        return mm_sum(r);
    }
    return mm_sum(r) + _mm512_reduce_add_pd(sum);
}

template<>
inline float
veryfasttree::AVX512Operations<float>::vector_dot_product_rot(float f1[], float f2[], float fBy[], int64_t n) {
    __m512 sum1 = _mm512_setzero_ps();
    __m512 sum2 = _mm512_setzero_ps();
    int64_t m = n - (n % 16);
    for (int64_t i = 0; i < m; i += 16) {
        __m512 a1, a2, aBy;
        a1 = _mm512_load_ps(f1 + i);
        a2 = _mm512_load_ps(f2 + i);
        aBy = _mm512_load_ps(fBy + i);
        sum1 = _mm512_add_ps(_mm512_mul_ps(a1, aBy), sum1);
        sum2 = _mm512_add_ps(_mm512_mul_ps(a2, aBy), sum2);
    }
    __m128 a1, a2, aBy, r1, r2;
    a1 = _mm_load_ps(f1 + m);
    a2 = _mm_load_ps(f2 + m);
    aBy = _mm_load_ps(fBy + m);
    r1 = _mm_mul_ps(a1, aBy);
    r2 = _mm_mul_ps(a2, aBy);
    if (n == 4) {
        return mm_sum(r1) * mm_sum(r2);
    }
    return (_mm512_reduce_add_ps(sum1) + mm_sum(r1)) * (_mm512_reduce_add_ps(sum2) + mm_sum(r2));
}

template<>
inline double veryfasttree::AVX512Operations<double>::
vector_dot_product_rot(double f1[], double f2[], double fBy[], int64_t n) {
    __m512d sum1 = _mm512_setzero_pd();
    __m512d sum2 = _mm512_setzero_pd();
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m512d a1, a2, aBy;
        a1 = _mm512_load_pd(f1 + i);
        a2 = _mm512_load_pd(f2 + i);
        aBy = _mm512_load_pd(fBy + i);
        sum1 = _mm512_add_pd(_mm512_mul_pd(a1, aBy), sum1);
        sum2 = _mm512_add_pd(_mm512_mul_pd(a2, aBy), sum2);
    }
    __m256d a1, a2, aBy, r1, r2;
    a1 = _mm256_load_pd(f1 + m);
    a2 = _mm256_load_pd(f2 + m);
    aBy = _mm256_load_pd(fBy + m);
    r1 = _mm256_mul_pd(a1, aBy);
    r2 = _mm256_mul_pd(a2, aBy);
    if (n == 4) {
        return mm_sum(r1) * mm_sum(r2);
    }
    return (_mm512_reduce_add_pd(sum1) + mm_sum(r1)) * (_mm512_reduce_add_pd(sum2) + mm_sum(r2));
}

template<>
inline void veryfasttree::AVX512Operations<float>::vector_add(float fTot[], float fAdd[], int64_t n) {
    __m512 a, b;
    for (int64_t i = 0; i < n; i += 16) {
        a = _mm512_load_ps(fTot + i);
        b = _mm512_load_ps(fAdd + i);
        _mm512_store_ps(fTot + i, _mm512_add_ps(a, b));
    }
}

template<>
inline void veryfasttree::AVX512Operations<double>::vector_add(double fTot[], double fAdd[], int64_t n) {
    __m512d a, b;
    for (int64_t i = 0; i < n; i += 8) {
        a = _mm512_load_pd(fTot + i);
        b = _mm512_load_pd(fAdd + i);
        _mm512_store_pd(fTot + i, _mm512_add_pd(a, b));
    }
}

template<>
inline float veryfasttree::AVX512Operations<float>::vector_sum(float f1[], int64_t n) {
    __m512 sum = _mm512_setzero_ps();
    int64_t m = n - (n % 16);
    for (int64_t i = 0; i < m; i += 16) {
        __m512 a = _mm512_load_ps(f1 + i);
        sum = _mm512_add_ps(a, sum);
    }
    __m128 a = _mm_load_ps(f1 + m);
    if (n == 4) {
        return mm_sum(a);
    }
    return _mm512_reduce_add_ps(sum) + mm_sum(a);
}

template<>
inline double veryfasttree::AVX512Operations<double>::vector_sum(double f1[], int64_t n) {
    __m512d sum = _mm512_setzero_pd();
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m512d a = _mm512_load_pd(f1 + i);
        sum = _mm512_add_pd(a, sum);
    }
    __m256d a = _mm256_load_pd(f1 + m);
    if (n == 4) {
        return mm_sum(a);
    }
    return _mm512_reduce_add_pd(sum) + mm_sum(a);
}

template<>
inline void
veryfasttree::AVX512Operations<float>::vector_multiply_by(float f[], float fBy, int64_t n, float fOut[]) {
    __m512 c = _mm512_set1_ps(fBy);
    int64_t m = n - (n % 16);
    for (int64_t i = 0; i < m; i += 16) {
        __m512 a, b;
        a = _mm512_load_ps(f + i);
        b = _mm512_mul_ps(a, c);
        _mm512_store_ps(fOut + i, b);
    }
    __m128 a, b;
    a = _mm_load_ps(f + m);
    b = _mm_mul_ps(a, _mm512_castps512_ps128(c));
    _mm_store_ps(fOut + m, b);
}

template<>
inline void
veryfasttree::AVX512Operations<double>::vector_multiply_by(double f[], double fBy, int64_t n, double fOut[]) {
    __m512d c = _mm512_set1_pd(fBy);
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m512d a, b;
        a = _mm512_load_pd(f + i);
        b = _mm512_mul_pd(a, c);
        _mm512_store_pd(fOut + i, b);
    }
    __m256d a, b;
    a = _mm256_load_pd(f + m);
    b = _mm256_mul_pd(a, _mm512_castpd512_pd256(c));
    _mm256_store_pd(fOut + m, b);
}

template<>
inline void
veryfasttree::AVX512Operations<float>::vector_add_mult(float fTot[], float fAdd[], float weight, int64_t n) {
    __m512 w = _mm512_set1_ps(weight);
    int64_t m = n - (n % 16);
    for (int64_t i = 0; i < m; i += 16) {
        __m512 tot, add;
        tot = _mm512_load_ps(fTot + i);
        add = _mm512_load_ps(fAdd + i);
        _mm512_store_ps(fTot + i, _mm512_add_ps(tot, _mm512_mul_ps(add, w)));
    }
    __m128 tot, add;
    tot = _mm_load_ps(fTot + m);
    add = _mm_load_ps(fAdd + m);
    _mm_store_ps(fTot + m, _mm_add_ps(tot, _mm_mul_ps(add, _mm512_castps512_ps128(w))));
}

template<>
inline void
veryfasttree::AVX512Operations<double>::vector_add_mult(double fTot[], double fAdd[], double weight, int64_t n) {
    __m512d w = _mm512_set1_pd(weight);
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m512d tot, add;
        tot = _mm512_load_pd(fTot + i);
        add = _mm512_load_pd(fAdd + i);
        _mm512_store_pd(fTot + i, _mm512_add_pd(tot, _mm512_mul_pd(add, w)));
    }
    __m256d tot, add;
    tot = _mm256_load_pd(fTot + m);
    add = _mm256_load_pd(fAdd + m);
    _mm256_store_pd(fTot + m, _mm256_add_pd(tot, _mm256_mul_pd(add, _mm512_castpd512_pd256(w))));
}

template<>
template<int row>
inline void veryfasttree::AVX512Operations<float>::matrix_by_vector4(float mat[][row], float vec[], float out[]) {
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
template<int row>
inline void
veryfasttree::AVX512Operations<double>::matrix_by_vector4(double mat[][row], double vec[], double out[]) {
    __m256d o = _mm256_setzero_pd();
    for (int64_t j = 0; j < 4; j++) {
        __m256d m = _mm256_load_pd(&mat[j][0]);
        __m256d vj = _mm256_broadcast_sd(&vec[j]);
        o = _mm256_add_pd(o, _mm256_mul_pd(vj, m));
    }
    _mm256_store_pd(out, o);
}

template<>
inline void veryfasttree::AVX512Operations<float>::fastexp(float fTot[], int64_t n, int lvl) {
    if (lvl == 0) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = (float) std::exp((double) fTot[k]);
        }
    } else if (lvl == 1) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = std::exp(fTot[k]);
        }
    } else if (lvl == 2) {
        alignas(ALIGNMENT) double arr[8];
        int64_t m = n - (n % 8);
        for (int64_t k = 0; k < m; k += 8) {
            __m512d nums = _mm512_set_pd(fTot[k + 7], fTot[k + 6], fTot[k + 5], fTot[k + 4],
                                         fTot[k + 3], fTot[k + 2], fTot[k + 1], fTot[k]);
            __m512d res = fastexpImpl(nums);
            _mm512_store_pd(arr, res);
            fTot[k] = (float) arr[0];
            fTot[k + 1] = (float) arr[1];
            fTot[k + 2] = (float) arr[2];
            fTot[k + 3] = (float) arr[3];
            fTot[k + 4] = (float) arr[4];
            fTot[k + 5] = (float) arr[5];
            fTot[k + 6] = (float) arr[6];
            fTot[k + 7] = (float) arr[7];
        }
        __m512d nums = _mm512_set_pd(fTot[m + 3], fTot[m + 2], fTot[m + 1], fTot[m],
                                     fTot[m + 3], fTot[m + 2], fTot[m + 1], fTot[m]);
        __m512d res = fastexpImpl(nums);
        _mm512_store_pd(arr, res);
        fTot[m] = (float) arr[0];
        fTot[m + 1] = (float) arr[1];
        fTot[m + 2] = (float) arr[2];
        fTot[m + 3] = (float) arr[3];
    } else {
        int64_t m = n - (n % 16);
        for (int64_t k = 0; k < m; k += 16) {
            __m512 nums = _mm512_load_ps(&fTot[k]);
            __m512 res = fastexpImpl(nums);
            _mm512_store_ps(&fTot[k], res);
        }
        __m128 nums = _mm_load_ps(&fTot[m]);
        __m128 res = fastexpImpl(nums);
        _mm_store_ps(&fTot[m], res);
    }
}

template<>
inline void veryfasttree::AVX512Operations<double>::fastexp(double fTot[], int64_t n, int lvl) {
    if (lvl == 0) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = std::exp(fTot[k]);
        }
    } else if (lvl == 1) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = (double) std::exp((float) fTot[k]);
        }
    } else if (lvl == 2) {
        int64_t m = n - (n % 16);
        for (int64_t k = 0; k < m; k += 8) {
            __m512d nums = _mm512_load_pd(&fTot[k]);
            __m512d res = fastexpImpl(nums);
            _mm512_store_pd(&fTot[k], res);
        }
        __m256d nums = _mm256_load_pd(&fTot[m]);
        __m256d res = fastexpImpl(nums);
        _mm256_store_pd(&fTot[m], res);
    } else {
        int64_t m = n - (n % 16);
        alignas(ALIGNMENT) float arr[16];
        for (int64_t k = 0; k < m; k += 16) {
            __m512 nums = _mm512_set_ps((float) fTot[k + 15], (float) fTot[k + 14], (float) fTot[k + 13],
                                        (float) fTot[k + 12], (float) fTot[k + 11], (float) fTot[k + 10],
                                        (float) fTot[k + 9], (float) fTot[k + 8], (float) fTot[k + 7],
                                        (float) fTot[k + 6], (float) fTot[k + 5], (float) fTot[k + 4],
                                        (float) fTot[k + 3], (float) fTot[k + 2], (float) fTot[k + 1], (float) fTot[k]);
            __m512 res = fastexpImpl(nums);
            _mm512_store_ps(arr, res);
            fTot[k] = (float) arr[0];
            fTot[k + 1] = (float) arr[1];
            fTot[k + 2] = (float) arr[2];
            fTot[k + 3] = (float) arr[3];
            fTot[k + 4] = (float) arr[4];
            fTot[k + 5] = (float) arr[5];
            fTot[k + 6] = (float) arr[6];
            fTot[k + 7] = (float) arr[7];
            fTot[k + 8] = (float) arr[8];
            fTot[k + 9] = (float) arr[9];
            fTot[k + 10] = (float) arr[10];
            fTot[k + 11] = (float) arr[11];
            fTot[k + 12] = (float) arr[12];
            fTot[k + 13] = (float) arr[13];
            fTot[k + 14] = (float) arr[14];
            fTot[k + 15] = (float) arr[15];
        }
        __m128 nums = _mm_set_ps((float) fTot[m + 3], (float) fTot[m + 2], (float) fTot[m + 1], (float) fTot[m]);
        __m128 res = fastexpImpl(nums);
        _mm_store_ps(arr, res);
        fTot[m] = (float) arr[0];
        fTot[m + 1] = (float) arr[1];
        fTot[m + 2] = (float) arr[2];
        fTot[m + 3] = (float) arr[3];
    }
}

template<typename Precision>
inline __m512 veryfasttree::AVX512Operations<Precision>::fastexpImpl(__m512 x) {
    __m512i m, u;
    __m512 xx, px, qx;

    px = _mm512_add_ps(_mm512_mul_ps(_mm512_set1_ps(1.4426950408889634073599f), x), _mm512_set1_ps(0.5));
    px = _mm512_floor_ps(px);

    /* m = (__m128i) px */
    m = _mm512_cvttps_epi32(px);

    x = _mm512_sub_ps(x, _mm512_mul_ps(px, _mm512_set1_ps(6.93145751953125E-1f)));
    x = _mm512_sub_ps(x, _mm512_mul_ps(px, _mm512_set1_ps(1.42860682030941723212E-6f)));
    xx = _mm512_mul_ps(x, x);

    /* px = x * P(x**2). */
    px = _mm512_set1_ps(1.26177193074810590878E-4f);
    px = _mm512_mul_ps(px, xx);
    px = _mm512_add_ps(px, _mm512_set1_ps(3.02994407707441961300E-2f));
    px = _mm512_mul_ps(px, xx);
    px = _mm512_add_ps(px, _mm512_set1_ps(9.99999999999999999910E-1f));
    px = _mm512_mul_ps(px, x);

    /* Evaluate Q(x**2). */
    qx = _mm512_set1_ps(3.00198505138664455042E-6f);
    qx = _mm512_mul_ps(qx, xx);
    qx = _mm512_add_ps(qx, _mm512_set1_ps(2.52448340349684104192E-3f));
    qx = _mm512_mul_ps(qx, xx);
    qx = _mm512_add_ps(qx, _mm512_set1_ps(2.27265548208155028766E-1f));
    qx = _mm512_mul_ps(qx, xx);
    qx = _mm512_add_ps(qx, _mm512_set1_ps(2.00000000000000000009E0f));

    /* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) ) */
    x = _mm512_div_ps(px, _mm512_sub_ps(qx, px));
    x = _mm512_add_ps(_mm512_set1_ps(1.0), _mm512_mul_ps(_mm512_set1_ps(2.0), x));

    /* Build 2^n in float. */
    m = _mm512_add_epi32(m, _mm512_set1_epi32(127));
    u = _mm512_slli_epi32(m, 23);

    return _mm512_mul_ps(x, _mm512_castsi512_ps(u));
}

template<typename Precision>
inline __m512d veryfasttree::AVX512Operations<Precision>::fastexpImpl(__m512d x) {
    __m512i u;
    __m512d xx, px, qx;

    px = _mm512_add_pd(_mm512_mul_pd(_mm512_set1_pd(1.4426950408889634073599), x), _mm512_set1_pd(0.5));
    px = _mm512_floor_pd(px);

    /* m = (__m512i) px */
    __m512i m = _mm512_sub_epi64(
            _mm512_castpd_si512(_mm512_add_pd(px, _mm512_set1_pd(0x0018000000000000))),
            _mm512_castpd_si512(_mm512_set1_pd(0x0018000000000000))
    );

    /* Build 2^n in double. */
    m = _mm512_add_epi64(m, _mm512_set1_epi64(1023));
    u = _mm512_slli_epi64(m, 52);

    x = _mm512_sub_pd(x, _mm512_mul_pd(px, _mm512_set1_pd(6.93145751953125E-1)));
    x = _mm512_sub_pd(x, _mm512_mul_pd(px, _mm512_set1_pd(1.42860682030941723212E-6)));
    xx = _mm512_mul_pd(x, x);

    /* px = x * P(x**2). */
    px = _mm512_set1_pd(1.26177193074810590878E-4);
    px = _mm512_mul_pd(px, xx);
    px = _mm512_add_pd(px, _mm512_set1_pd(3.02994407707441961300E-2));
    px = _mm512_mul_pd(px, xx);
    px = _mm512_add_pd(px, _mm512_set1_pd(9.99999999999999999910E-1));
    px = _mm512_mul_pd(px, x);

    /* Evaluate Q(x**2). */
    qx = _mm512_set1_pd(3.00198505138664455042E-6);
    qx = _mm512_mul_pd(qx, xx);
    qx = _mm512_add_pd(qx, _mm512_set1_pd(2.52448340349684104192E-3));
    qx = _mm512_mul_pd(qx, xx);
    qx = _mm512_add_pd(qx, _mm512_set1_pd(2.27265548208155028766E-1));
    qx = _mm512_mul_pd(qx, xx);
    qx = _mm512_add_pd(qx, _mm512_set1_pd(2.00000000000000000009E0));

    /* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) ) */
    x = _mm512_div_pd(px, _mm512_sub_pd(qx, px));
    x = _mm512_add_pd(_mm512_set1_pd(1.0), _mm512_mul_pd(_mm512_set1_pd(2.0), x));

    return _mm512_mul_pd(x, _mm512_castsi512_pd(u));
}

/* Copy from SSE256Operations.tcc */
template<typename Precision>
inline __m256d veryfasttree::AVX512Operations<Precision>::fastexpImpl(__m256d x) {
    __m256i u;
    __m256d xx, px, qx;

    px = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(1.4426950408889634073599), x), _mm256_set1_pd(0.5));
    px = _mm256_floor_pd(px);

    #ifdef __AVX2__
    /* m = (__m256i) px */
    __m256i m = _mm256_sub_epi64(
            _mm256_castpd_si256(_mm256_add_pd(px, _mm256_set1_pd(0x0018000000000000))),
            _mm256_castpd_si256(_mm256_set1_pd(0x0018000000000000))
    );

    /* Build 2^n in double. */
    m = _mm256_add_epi64(m, _mm256_set1_epi64x(1023));
    u = _mm256_slli_epi64(m, 52);
    #else
    /* m = (__m256i) px */
    __m128i lm = _mm_sub_epi64(
            _mm_castpd_si128(_mm_add_pd(_mm256_castpd256_pd128(px), _mm_set1_pd(0x0018000000000000))),
            _mm_castpd_si128(_mm_set1_pd(0x0018000000000000))
    );
    __m128i hm = _mm_sub_epi64(
            _mm_castpd_si128(_mm_add_pd(_mm256_extractf128_pd(px, 1), _mm_set1_pd(0x0018000000000000))),
            _mm_castpd_si128(_mm_set1_pd(0x0018000000000000))
    );

    /* Build 2^n in double. */
    lm = _mm_add_epi64(lm, _mm_set1_epi64x(1023));
    hm = _mm_add_epi64(hm, _mm_set1_epi64x(1023));
    __m128i lu = _mm_slli_epi64(lm, 52);
    __m128i hu = _mm_slli_epi64(hm, 52);
    u = _mm256_castsi128_si256(lu);
    u = _mm256_insertf128_si256(u, hu, 1);
    #endif

    x = _mm256_sub_pd(x, _mm256_mul_pd(px, _mm256_set1_pd(6.93145751953125E-1)));
    x = _mm256_sub_pd(x, _mm256_mul_pd(px, _mm256_set1_pd(1.42860682030941723212E-6)));
    xx = _mm256_mul_pd(x, x);

    /* px = x * P(x**2). */
    px = _mm256_set1_pd(1.26177193074810590878E-4);
    px = _mm256_mul_pd(px, xx);
    px = _mm256_add_pd(px, _mm256_set1_pd(3.02994407707441961300E-2));
    px = _mm256_mul_pd(px, xx);
    px = _mm256_add_pd(px, _mm256_set1_pd(9.99999999999999999910E-1));
    px = _mm256_mul_pd(px, x);

    /* Evaluate Q(x**2). */
    qx = _mm256_set1_pd(3.00198505138664455042E-6);
    qx = _mm256_mul_pd(qx, xx);
    qx = _mm256_add_pd(qx, _mm256_set1_pd(2.52448340349684104192E-3));
    qx = _mm256_mul_pd(qx, xx);
    qx = _mm256_add_pd(qx, _mm256_set1_pd(2.27265548208155028766E-1));
    qx = _mm256_mul_pd(qx, xx);
    qx = _mm256_add_pd(qx, _mm256_set1_pd(2.00000000000000000009E0));

    /* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) ) */
    x = _mm256_div_pd(px, _mm256_sub_pd(qx, px));
    x = _mm256_add_pd(_mm256_set1_pd(1.0), _mm256_mul_pd(_mm256_set1_pd(2.0), x));

    return _mm256_mul_pd(x, _mm256_castsi256_pd(u));
}

/* Copy from SSE128Operations.tcc */
template<typename Precision>
inline __m128 veryfasttree::AVX512Operations<Precision>::fastexpImpl(__m128 x) {
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

#undef getf
#undef getd
