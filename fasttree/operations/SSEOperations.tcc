
#ifndef FASTTREE_SSEOPERATIONS_TCC
#define FASTTREE_SSEOPERATIONS_TCC

#include "SSEOperations.h"

template<>
void fasttree::SSEOperations<float>::
vector_multiply(float  f1[], float  f2[], int64_t n, float  fOut[]) {
    for (int64_t i = 0; i < n; i += 4) {
        __m128 a, b, c;
        a = _mm_load_ps(f1 + i);
        b = _mm_load_ps(f2 + i);
        c = _mm_mul_ps(a, b);
        _mm_store_ps(fOut + i, c);
    }
}

template<>
void fasttree::SSEOperations<double>::
vector_multiply(double  f1[], double  f2[], int64_t n, double  fOut[]) {
    for (int64_t i = 0; i < n; i += 2) {
        __m128d a, b, c;
        a = _mm_load_pd(f1 + i);
        b = _mm_load_pd(f2 + i);
        c = _mm_mul_pd(a, b);
        _mm_store_pd(fOut + i, c);
    }
}

template<>
float fasttree::SSEOperations<float>
::vector_multiply_sum(float  f1[], float  f2[], int64_t n) {
    if (n == 4) {
        return (f1[0] * f2[0] + f1[1] * f2[1] + f1[2] * f2[2] + f1[3] * f2[3]);
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
double fasttree::SSEOperations<double>
::vector_multiply_sum(double  f1[], double  f2[], int64_t n) {
    if (n == 4) {
        return (f1[0] * f2[0] + f1[1] * f2[1] + f1[2] * f2[2] + f1[3] * f2[3]);
    }
    __m128d sum = _mm_setzero_pd();
    for (int64_t i = 0; i < n; i += 4) {
        __m128d a, b, c;
        a = _mm_load_pd(f1 + i);
        b = _mm_load_pd(f2 + i);
        c = _mm_mul_pd(a, b);
        sum = _mm_add_pd(c, sum);
    }
    return mm_sum(sum);
}

template<>
void fasttree::SSEOperations<float>::
vector_add_mult(float  fTot[], float  fAdd[], float weight, int64_t n) {
    __m128 w = _mm_set1_ps(weight);
    for (int64_t i = 0; i < n; i += 4) {
        __m128 tot, add;
        tot = _mm_load_ps(fTot + i);
        add = _mm_load_ps(fAdd + i);
        _mm_store_ps(fTot + i, _mm_add_ps(tot, _mm_mul_ps(add, w)));
    }
}

template<>
void fasttree::SSEOperations<double>::
vector_add_mult(double  fTot[], double  fAdd[], double weight, int64_t n) {
    __m128d w = _mm_set1_pd(weight);
    for (int64_t i = 0; i < n; i += 2) {
        __m128d tot, add;
        tot = _mm_load_pd(fTot + i);
        add = _mm_load_pd(fAdd + i);
        _mm_store_pd(fTot + i, _mm_add_pd(tot, _mm_mul_pd(add, w)));
    }
}


template<>
float fasttree::SSEOperations<float>::mm_sum(register __m128 sum) {
    /* stupider but faster */
    alignas(ALIGNMENT) float f[4];
    _mm_store_ps(f, sum);
    return (f[0] + f[1] + f[2] + f[3]);
}

template<>
double fasttree::SSEOperations<double>::mm_sum(register __m128 sum) {
    /* stupider but faster */
    alignas(ALIGNMENT) double f[4];
    ALIGNED;
    _mm_store_pd(f, sum);
    return (f[0] + f[1] + f[2] + f[3]);
}

#endif

