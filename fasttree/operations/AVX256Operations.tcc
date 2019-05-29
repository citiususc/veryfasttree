
#ifndef FASTTREE_AVX256OPERATIONS_TCC
#define FASTTREE_AVX256OPERATIONS_TCC

#include "AVX256Operations.tcc"

template<>
template<>
inline double fasttree::AVX256Operations<double>::mm_sum(register __m256d sum) {
    alignas(ALIGNMENT) double f[4];
    _mm256_store_pd(f, sum);
    return (f[0] + f[1] + f[2] + f[3]);
}

template<>
template<>
inline float fasttree::AVX256Operations<float>::mm_sum(register __m256 sum) {
    __m256 sum4 = _mm256_hadd_ps(sum, sum);
    alignas(ALIGNMENT) float f[8];
    _mm256_store_ps(f, sum4);
    return (f[0] + f[1] + f[4] + f[5]);
}

template<>
void fasttree::AVX256Operations<double>::
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
void fasttree::AVX256Operations<float>::
vector_multiply(float f1[], float f2[], int64_t n, float fOut[]) {
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 a, b, c;
        a = _mm256_load_ps(f1 + i);
        b = _mm256_load_ps(f2 + i);
        c = _mm256_mul_ps(a, b);
        _mm256_store_ps(fOut + i, c);
    }
    if (n != m) {
        __m128 a, b, c;
        a = _mm_load_ps(f1 + m);
        b = _mm_load_ps(f2 + m);
        c = _mm_mul_ps(a, b);
        _mm_store_ps(fOut + m, c);
    }
}

template<>
double fasttree::AVX256Operations<double>
::vector_multiply_sum(double f1[], double f2[], int64_t n) {
    if (n == 4) {
        return f1[0] * f2[0] + f1[1] * f2[1] + f1[2] * f2[2] + f1[3] * f2[3];
    }
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
float fasttree::AVX256Operations<float>
::vector_multiply_sum(float f1[], float f2[], int64_t n) {
    if (n == 4) {
        return f1[0] * f2[0] + f1[1] * f2[1] + f1[2] * f2[2] + f1[3] * f2[3];
    }
    __m256 sum = _mm256_setzero_ps();
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 a, b, c;
        a = _mm256_load_ps(f1 + i);
        b = _mm256_load_ps(f2 + i);
        c = _mm256_mul_ps(a, b);
        sum = _mm256_add_ps(c, sum);
    }
    if (n != m) {
        __m128 a, b, c;
        a = _mm_load_ps(f1 + m);
        b = _mm_load_ps(f2 + m);
        c = _mm_mul_ps(a, b);
        sum = _mm256_add_ps(_mm256_castps128_ps256(c), sum);
    }
    return mm_sum(sum);
}

template<>
double fasttree::AVX256Operations<double>::
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
float fasttree::AVX256Operations<float>::
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
    if (n != m) {
        __m128 a1, a2, a3, r;
        a1 = _mm_load_ps(f1 + m);
        a2 = _mm_load_ps(f2 + m);
        a3 = _mm_load_ps(f3 + m);
        r = _mm_mul_ps(_mm_mul_ps(a1, a2), a3);
        sum = _mm256_add_ps(_mm256_castps128_ps256(r), sum);
    }
    return mm_sum(sum);
}

template<>
double fasttree::AVX256Operations<double>::
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
    return (mm_sum(sum1) * mm_sum(sum2));
}

template<>
float fasttree::AVX256Operations<float>::
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
    if (n != m) {
        __m128 a1, a2, aBy, r1, r2;
        a1 = _mm_load_ps(f1 + m);
        a2 = _mm_load_ps(f2 + m);
        aBy = _mm_load_ps(fBy + m);
        r1 = _mm_mul_ps(a1, aBy);
        r2 = _mm_mul_ps(a2, aBy);
        sum1 = _mm256_add_ps(_mm256_castps128_ps256(r1), sum1);
        sum2 = _mm256_add_ps(_mm256_castps128_ps256(r2), sum2);
    }
    return mm_sum(sum1) * mm_sum(sum2);
}

template<>
double fasttree::AVX256Operations<double>::
vector_sum(double f1[], int64_t n) {
    if (n == 4) {
        return f1[0] + f1[1] + f1[2] + f1[3];
    }
    __m256d sum = _mm256_setzero_pd();
    for (int64_t i = 0; i < n; i += 4) {
        __m256d a;
        a = _mm256_load_pd(f1 + i);
        sum = _mm256_add_pd(a, sum);
    }
    return (mm_sum(sum));
}

template<>
float fasttree::AVX256Operations<float>::
vector_sum(float f1[], int64_t n) {
    if (n == 4) {
        return f1[0] + f1[1] + f1[2] + f1[3];
    }
    __m256 sum = _mm256_setzero_ps();
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 a;
        a = _mm256_load_ps(f1 + i);
        sum = _mm256_add_ps(a, sum);
    }
    if (n != m) {
        __m128 a = _mm_load_ps(f1 + m);
        sum = _mm256_add_ps(_mm256_castps128_ps256(a), sum);
    }
    return mm_sum(sum);
}

template<>
void fasttree::AVX256Operations<double>::
vector_multiply_by(double f[], double fBy, int64_t n) {
    __m256d c = _mm256_set1_pd(fBy);
    for (int64_t i = 0; i < n; i += 4) {
        __m256d a, b;
        a = _mm256_load_pd(f + i);
        b = _mm256_mul_pd(a, c);
        _mm256_store_pd(f + i, b);
    }
}

template<>
void fasttree::AVX256Operations<float>::
vector_multiply_by(float f[], float fBy, int64_t n) {
    __m256 c = _mm256_set1_ps(fBy);
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 a, b;
        a = _mm256_load_ps(f + i);
        b = _mm256_mul_ps(a, c);
        _mm256_store_ps(f + i, b);
    }
    if (n != m) {
        __m128 a, b;
        a = _mm_load_ps(f + m);
        b = _mm_mul_ps(a, _mm256_castps256_ps128(c));
        _mm_store_ps(f + m, b);
    }
}

template<>
void fasttree::AVX256Operations<double>::
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
void fasttree::AVX256Operations<float>::
vector_add_mult(float fTot[], float fAdd[], float weight, int64_t n) {
    __m256 w = _mm256_set1_ps(weight);
    int64_t m = n - (n % 8);
    for (int64_t i = 0; i < m; i += 8) {
        __m256 tot, add;
        tot = _mm256_load_ps(fTot + i);
        add = _mm256_load_ps(fAdd + i);
        _mm256_store_ps(fTot + i, _mm256_add_ps(tot, _mm256_mul_ps(add, w)));
    }
    if (n != m) {
        __m128 tot, add;
        tot = _mm_load_ps(fTot + m);
        add = _mm_load_ps(fAdd + m);
        _mm_store_ps(fTot + m, _mm_add_ps(tot, _mm_mul_ps(add, _mm256_castps256_ps128(w))));
    }
}

template<>
void fasttree::AVX256Operations<double>::
matrixt_by_vector4(double mat[4][4], double vec[4], double out[4]) {
    __m256d o = _mm256_setzero_pd();
    /* result is a sum of vectors: sum(k) v[k] * mat[k][] */
    for (int64_t j = 0; j < 4; j++) {
        __m256d m = _mm256_load_pd(&mat[j][0]);
        __m256d vj = _mm256_broadcast_sd(&vec[j]);    /* is it faster to shuffle v? */
        o = _mm256_add_pd(o, _mm256_mul_pd(vj, m));
    }
    _mm256_store_pd(out, o);
}

template<>
void fasttree::AVX256Operations<float>::
matrixt_by_vector4(float mat[4][4], float vec[4], float out[4]) {
    __m256 vj1 = _mm256_set_ps(vec[0], vec[0], vec[0], vec[0], vec[1], vec[1], vec[1], vec[1]);
    __m256 vj2 = _mm256_set_ps(vec[2], vec[2], vec[2], vec[2], vec[3], vec[3], vec[3], vec[3]);
    __m256 m1 = _mm256_load_ps(&mat[0][0]);
    __m256 m2 = _mm256_load_ps(&mat[2][0]);
    __m256 o1 = _mm256_mul_ps(vj1, m1);
    __m256 o2 = _mm256_mul_ps(vj2, m2);
    __m256 s1 = _mm256_add_ps(o1, o2);
    __m128 s2 = _mm_add_ps(_mm256_extractf128_ps(s1, 0), _mm256_extractf128_ps(s1, 1));
    _mm_store1_ps(out, s2);
}

#endif
