
#ifndef FASTTREE_SSE3OPERATIONS_TCC
#define FASTTREE_SSE3OPERATIONS_TCC

#include "SSE3Operations.h"

template<>
template<>
inline float fasttree::SSE3Operations<float>::mm_sum(register __m128 sum) {
    /* stupider but faster */
    alignas(ALIGNMENT) float f[4];
    _mm_store_ps(f, sum);
    return (f[0] + f[1] + f[2] + f[3]);
}

template<>
template<>
inline double fasttree::SSE3Operations<double>::mm_sum(register __m128d sum) {
    /* stupider but faster */
    alignas(ALIGNMENT) double f[4];;
    _mm_store_pd(f, sum);
    return (f[0] + f[1] + f[2] + f[3]);
}

template<>
void fasttree::SSE3Operations<float>::
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
void fasttree::SSE3Operations<double>::
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
float fasttree::SSE3Operations<float>
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
double fasttree::SSE3Operations<double>
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
float fasttree::SSE3Operations<float>::
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
double fasttree::SSE3Operations<double>::
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
float fasttree::SSE3Operations<float>::
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
double fasttree::SSE3Operations<double>::
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
float fasttree::SSE3Operations<float>::
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
double fasttree::SSE3Operations<double>::
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
void fasttree::SSE3Operations<float>::
vector_multiply_by(float f[], float fBy, int64_t n) {
    __m128 c = _mm_set1_ps(fBy);
    for (int64_t i = 0; i < n; i += 4) {
        __m128 a, b;
        a = _mm_load_ps(f + i);
        b = _mm_mul_ps(a, c);
        _mm_store_ps(f + i, b);
    }
}

template<>
void fasttree::SSE3Operations<double>::
vector_multiply_by(double f[], double fBy, int64_t n) {
    __m128d c = _mm_set1_pd(fBy);
    for (int64_t i = 0; i < n; i += 2) {
        __m128d a, b;
        a = _mm_load_pd(f + i);
        b = _mm_mul_pd(a, c);
        _mm_store_pd(f + i, b);
    }
}

template<>
void fasttree::SSE3Operations<float>::
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
void fasttree::SSE3Operations<double>::
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
void fasttree::SSE3Operations<float>::
matrixt_by_vector4(numeric_t mat[4][4], numeric_t vec[4], numeric_t out[4]) {
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
void fasttree::SSE3Operations<double>::
matrixt_by_vector4(numeric_t mat[4][4], numeric_t vec[4], numeric_t out[4]) {
    __m128d o1 = _mm_setzero_pd();
    __m128d o2 = _mm_setzero_pd();
    /* result is a sum of vectors: sum(k) v[k] * mat[k][] */
    for (int64_t j = 0; j < 4; j++) {
        __m128d vj = _mm_load1_pd(&vec[j]); /* is it faster to shuffle v? */
        __m128d m1 = _mm_load_pd(&mat[j][0]);
        __m128d m2= _mm_load_pd(&mat[j][2]);
        o1 = _mm_add_pd(o1, _mm_mul_pd(vj, m1));
        o2 = _mm_add_pd(o2, _mm_mul_pd(vj, m2));
    }
    _mm_store_pd(out, o1);
    _mm_store_pd(&out[2], o2);
}

#endif

