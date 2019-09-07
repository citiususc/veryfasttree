
#ifndef FASTTREE_AVX256OPERATIONS_TCC
#define FASTTREE_AVX256OPERATIONS_TCC

#include "AVX256Operations.tcc"

template<>
inline float fasttree::AVX256Operations<float>::mm_sum(__m128 sum) {
    sum = _mm_hadd_ps(sum, sum);
    return sum[0] + sum[1];
}

template<>
inline double fasttree::AVX256Operations<double>::mm_sum(__m256d sum) {
    sum = _mm256_hadd_pd(sum, sum);
    return sum[0] + sum[2];
}

template<>
inline float fasttree::AVX256Operations<float>::mm_sum(__m256 sum1, __m128 sum2) {
    sum2 = _mm_add_ps(*(__m128 *) &sum1, sum2);
    sum2 = _mm_add_ps(*(__m128 *) &sum1[4], sum2);
    sum2 = _mm_hadd_ps(sum2, sum2);
    return sum2[0] + sum2[1];
}

template<>
inline void fasttree::AVX256Operations<double>::
vector_multiply(double f1[], double f2[], int64_t n, double fOut[]) {
    for (int64_t i = 0; i < n; i += 4) {
        __m256d a, b, c;
        a = _mm256_load_pd(f1 + i);
        b = _mm256_load_pd(f2 + i);
        c = _mm256_mul_pd(a, b);
        _mm256_store_pd(fOut + i, c);
    }
}

template<>
inline void fasttree::AVX256Operations<float>::
vector_multiply(float f1[], float f2[], int64_t n, float fOut[]) {
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 a, b, c;
        a = _mm256_load_ps(f1 + i);
        b = _mm256_load_ps(f2 + i);
        c = _mm256_mul_ps(a, b);
        _mm256_store_ps(fOut + i, c);
    }
    __m128 a, b, c;
    a = _mm_load_ps(f1 + m);
    b = _mm_load_ps(f2 + m);
    c = _mm_mul_ps(a, b);
    _mm_store_ps(fOut + m, c);
}

template<>
inline double fasttree::AVX256Operations<double>
::vector_multiply_sum(double f1[], double f2[], int64_t n) {
    __m256d sum = _mm256_setzero_pd();
    for (int64_t i = 0; i < n; i += 4) {
        __m256d a, b, c;
        a = _mm256_load_pd(f1 + i);
        b = _mm256_load_pd(f2 + i);
        c = _mm256_mul_pd(a, b);
        sum = _mm256_add_pd(c, sum);
    }
    return mm_sum(sum);
}

template<>
inline float fasttree::AVX256Operations<float>
::vector_multiply_sum(float f1[], float f2[], int64_t n) {
    __m256 sum = _mm256_setzero_ps();
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 a, b, c;
        a = _mm256_load_ps(f1 + i);
        b = _mm256_load_ps(f2 + i);
        c = _mm256_mul_ps(a, b);
        sum = _mm256_add_ps(c, sum);
    }
    __m128 a, b, c;
    a = _mm_load_ps(f1 + m);
    b = _mm_load_ps(f2 + m);
    c = _mm_mul_ps(a, b);
    if (n == 4) {
        return mm_sum(c);
    }
    return mm_sum(sum, c);
}

template<>
inline double fasttree::AVX256Operations<double>::
vector_multiply3_sum(double f1[], double f2[], double f3[], int64_t n) {
    __m256d sum = _mm256_setzero_pd();
    for (int64_t i = 0; i < n; i += 4) {
        __m256d a1, a2, a3;
        a1 = _mm256_load_pd(f1 + i);
        a2 = _mm256_load_pd(f2 + i);
        a3 = _mm256_load_pd(f3 + i);
        sum = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(a1, a2), a3), sum);
    }
    return mm_sum(sum);
}

template<>
inline float fasttree::AVX256Operations<float>::
vector_multiply3_sum(float f1[], float f2[], float f3[], int64_t n) {
    __m256 sum = _mm256_setzero_ps();
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 a1, a2, a3;
        a1 = _mm256_load_ps(f1 + i);
        a2 = _mm256_load_ps(f2 + i);
        a3 = _mm256_load_ps(f3 + i);
        sum = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(a1, a2), a3), sum);
    }
    __m128 a1, a2, a3, r;
    a1 = _mm_load_ps(f1 + m);
    a2 = _mm_load_ps(f2 + m);
    a3 = _mm_load_ps(f3 + m);
    r = _mm_mul_ps(_mm_mul_ps(a1, a2), a3);
    if (n == 4) {
        return mm_sum(r);
    }
    return mm_sum(sum, r);
}

template<>
inline double fasttree::AVX256Operations<double>::
vector_dot_product_rot(double f1[], double f2[], double fBy[], int64_t n) {
    __m256d sum1 = _mm256_setzero_pd();
    __m256d sum2 = _mm256_setzero_pd();
    for (int64_t i = 0; i < n; i += 4) {
        __m256d a1, a2, aBy;
        a1 = _mm256_load_pd(f1 + i);
        a2 = _mm256_load_pd(f2 + i);
        aBy = _mm256_load_pd(fBy + i);
        sum1 = _mm256_add_pd(_mm256_mul_pd(a1, aBy), sum1);
        sum2 = _mm256_add_pd(_mm256_mul_pd(a2, aBy), sum2);
    }
    return mm_sum(sum1) * mm_sum(sum2);
}

template<>
inline float fasttree::AVX256Operations<float>::
vector_dot_product_rot(float f1[], float f2[], float fBy[], int64_t n) {
    __m256 sum1 = _mm256_setzero_ps();
    __m256 sum2 = _mm256_setzero_ps();
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 a1, a2, aBy;
        a1 = _mm256_load_ps(f1 + i);
        a2 = _mm256_load_ps(f2 + i);
        aBy = _mm256_load_ps(fBy + i);
        sum1 = _mm256_add_ps(_mm256_mul_ps(a1, aBy), sum1);
        sum2 = _mm256_add_ps(_mm256_mul_ps(a2, aBy), sum2);
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
    return mm_sum(sum1, r1) * mm_sum(sum2, r2);
}

template<>
inline double fasttree::AVX256Operations<double>::
vector_sum(double f1[], int64_t n) {
    __m256d sum = _mm256_setzero_pd();
    for (int64_t i = 0; i < n; i += 4) {
        __m256d a;
        a = _mm256_load_pd(f1 + i);
        sum = _mm256_add_pd(a, sum);
    }
    return mm_sum(sum);
}

template<>
inline float fasttree::AVX256Operations<float>::
vector_sum(float f1[], int64_t n) {
    __m256 sum = _mm256_setzero_ps();
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 a = _mm256_load_ps(f1 + i);
        sum = _mm256_add_ps(a, sum);
    }
    __m128 a = _mm_load_ps(f1 + m);
    if (n == 4) {
        return mm_sum(a);
    }
    return mm_sum(sum, a);
}

template<>
inline void fasttree::AVX256Operations<double>::
vector_multiply_by(double f[], double fBy, int64_t n, numeric_t fTot[]) {
    __m256d c = _mm256_set1_pd(fBy);
    for (int64_t i = 0; i < n; i += 4) {
        __m256d a, b;
        a = _mm256_load_pd(f + i);
        b = _mm256_mul_pd(a, c);
        _mm256_store_pd(fTot + i, b);
    }
}

template<>
inline void fasttree::AVX256Operations<float>::
vector_multiply_by(float f[], float fBy, int64_t n, numeric_t fTot[]) {
    __m256 c = _mm256_set1_ps(fBy);
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 a, b;
        a = _mm256_load_ps(f + i);
        b = _mm256_mul_ps(a, c);
        _mm256_store_ps(fTot + i, b);
    }
    __m128 a, b;
    a = _mm_load_ps(f + m);
    b = _mm_mul_ps(a, _mm256_castps256_ps128(c));
    _mm_store_ps(fTot + m, b);
}

template<>
inline void fasttree::AVX256Operations<double>::
vector_add_mult(double fTot[], double fAdd[], double weight, int64_t n) {
    __m256d w = _mm256_set1_pd(weight);
    for (int64_t i = 0; i < n; i += 4) {
        __m256d tot, add;
        tot = _mm256_load_pd(fTot + i);
        add = _mm256_load_pd(fAdd + i);
        _mm256_store_pd(fTot + i, _mm256_add_pd(tot, _mm256_mul_pd(add, w)));
    }
}

template<>
inline void fasttree::AVX256Operations<float>::
vector_add_mult(float fTot[], float fAdd[], float weight, int64_t n) {
    __m256 w = _mm256_set1_ps(weight);
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 tot, add;
        tot = _mm256_load_ps(fTot + i);
        add = _mm256_load_ps(fAdd + i);
        _mm256_store_ps(fTot + i, _mm256_add_ps(tot, _mm256_mul_ps(add, w)));
    }
    __m128 tot, add;
    tot = _mm_load_ps(fTot + m);
    add = _mm_load_ps(fAdd + m);
    _mm_store_ps(fTot + m, _mm_add_ps(tot, _mm_mul_ps(add, _mm256_castps256_ps128(w))));
}

template<>
template<int row>
inline void fasttree::AVX256Operations<double>::
matrixt_by_vector4(double mat[][row], double vec[], double out[]) {
    __m256d o = _mm256_setzero_pd();
    for (int64_t j = 0; j < 4; j++) {
        __m256d m = _mm256_load_pd(&mat[j][0]);
        __m256d vj = _mm256_broadcast_sd(&vec[j]);
        o = _mm256_add_pd(o, _mm256_mul_pd(vj, m));
    }
    _mm256_store_pd(out, o);
}

template<>
template<int row>
inline void fasttree::AVX256Operations<float>::
matrixt_by_vector4(float mat[][row], float vec[], float out[]) {
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
inline void fasttree::AVX256Operations<float>::fastexp(float fTot[], int64_t n, int lvl) {
    if (lvl == 0) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = (float) std::exp((double) fTot[k]);
        }
    } else if (lvl == 1) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = std::exp(fTot[k]);
        }
    } else if (lvl == 2) {
        alignas(ALIGNMENT) double arr[4];
        for (int64_t k = 0; k < n; k += 4) {
            __m256d nums = _mm256_set_pd(fTot[k + 3], fTot[k + 2], fTot[k + 1], fTot[k]);
            __m256d res = fastexpImpl(nums);
            _mm256_store_pd(arr, res);
            fTot[k] = (float) arr[0];
            fTot[k + 1] = (float) arr[1];
            fTot[k + 2] = (float) arr[2];
            fTot[k + 3] = (float) arr[3];
        }
    } else {
        int64_t m = n - (n % 8);
        for (int64_t k = 0; k < m; k += 8) {
            __m256 nums = _mm256_load_ps(&fTot[k]);
            __m256 res = fastexpImpl(nums);
            _mm256_store_ps(&fTot[k], res);
        }
        if (n != m) {
            __m128 nums = _mm_load_ps(&fTot[m]);
            __m128 res = fastexpImpl(nums);
            _mm_store_ps(&fTot[m], res);
        }
    }
}

template<>
inline void fasttree::AVX256Operations<double>::fastexp(double fTot[], int64_t n, int lvl) {
    if (lvl == 0) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = std::exp(fTot[k]);
        }
    } else if (lvl == 1) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = (double) std::exp((float) fTot[k]);
        }
    } else if (lvl == 2) {
        for (int64_t k = 0; k < n; k += 4) {
            __m256d nums = _mm256_load_pd(&fTot[k]);
            __m256d res = fastexpImpl(nums);
            _mm256_store_pd(&fTot[k], res);
        }
    } else {
        int64_t m = n - (n % 8);
        alignas(ALIGNMENT) float arr[8];
        for (int64_t k = 0; k < n; k += 8) {
            __m256 nums = _mm256_set_ps((float) fTot[k + 7], (float) fTot[k + 6], (float) fTot[k + 5],
                                        (float) fTot[k + 4],
                                        (float) fTot[k + 3], (float) fTot[k + 2], (float) fTot[k + 1], (float) fTot[k]);
            __m256 res = fastexpImpl(nums);
            _mm256_store_ps(arr, res);
            fTot[k] = (float) arr[0];
            fTot[k + 1] = (float) arr[1];
            fTot[k + 2] = (float) arr[2];
            fTot[k + 3] = (float) arr[3];
            fTot[k + 4] = (float) arr[4];
            fTot[k + 5] = (float) arr[5];
            fTot[k + 6] = (float) arr[6];
            fTot[k + 7] = (float) arr[7];
        }
        if (n != m) {
            __m128 nums = _mm_set_ps((float) fTot[m + 3], (float) fTot[m + 2], (float) fTot[m + 1], (float) fTot[m]);
            __m128 res = fastexpImpl(nums);
            _mm_store_ps(arr, res);
            fTot[m] = (float) arr[0];
            fTot[m + 1] = (float) arr[1];
            fTot[m + 2] = (float) arr[2];
            fTot[m + 3] = (float) arr[3];
        }
    }
}

template<typename Precision>
inline __m256 fasttree::AVX256Operations<Precision>::fastexpImpl(__m256 x) {
    __m256i m, u;
    __m256 xx, px, qx;

    px = _mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(1.4426950408889634073599f), x), _mm256_set1_ps(0.5));
    px = _mm256_floor_ps(px);

    /* m = (__m128i) px */
    m = _mm256_cvttps_epi32(px);

    x = _mm256_sub_ps(x, _mm256_mul_ps(px, _mm256_set1_ps(6.93145751953125E-1f)));
    x = _mm256_sub_ps(x, _mm256_mul_ps(px, _mm256_set1_ps(1.42860682030941723212E-6f)));
    xx = _mm256_mul_ps(x, x);

    /* px = x * P(x**2). */
    px = _mm256_set1_ps(1.26177193074810590878E-4f);
    px = _mm256_mul_ps(px, xx);
    px = _mm256_add_ps(px, _mm256_set1_ps(3.02994407707441961300E-2f));
    px = _mm256_mul_ps(px, xx);
    px = _mm256_add_ps(px, _mm256_set1_ps(9.99999999999999999910E-1f));
    px = _mm256_mul_ps(px, x);

    /* Evaluate Q(x**2). */
    qx = _mm256_set1_ps(3.00198505138664455042E-6f);
    qx = _mm256_mul_ps(qx, xx);
    qx = _mm256_add_ps(qx, _mm256_set1_ps(2.52448340349684104192E-3f));
    qx = _mm256_mul_ps(qx, xx);
    qx = _mm256_add_ps(qx, _mm256_set1_ps(2.27265548208155028766E-1f));
    qx = _mm256_mul_ps(qx, xx);
    qx = _mm256_add_ps(qx, _mm256_set1_ps(2.00000000000000000009E0f));

    /* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) ) */
    x = _mm256_div_ps(px, _mm256_sub_ps(qx, px));
    x = _mm256_add_ps(_mm256_set1_ps(1.0), _mm256_mul_ps(_mm256_set1_ps(2.0), x));

    /* Build 2^n in double. */
    #ifdef __AVX2__
    m = _mm256_add_epi32(m, _mm256_set1_epi32(127));
    u = _mm256_slli_epi32(m, 23);
    #else
    __m128i lm = _mm_add_epi32(_mm256_castsi256_si128(m), _mm_set1_epi32(127));
    __m128i hm = _mm_add_epi32(_mm256_extractf128_si256(m, 1), _mm_set1_epi32(127));
    __m128i lu = _mm_slli_epi32(lm, 23);
    __m128i hu = _mm_slli_epi32(hm, 23);
    u = _mm256_castsi128_si256(lu);
    u = _mm256_insertf128_si256(u, hu, 1);
    #endif

    return _mm256_mul_ps(x, _mm256_castsi256_ps(u));
}

template<typename Precision>
inline __m256d fasttree::AVX256Operations<Precision>::fastexpImpl(__m256d x) {
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
inline __m128 fasttree::AVX256Operations<Precision>::fastexpImpl(__m128 x) {
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

#endif
