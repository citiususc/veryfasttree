#ifndef FASTTREE_SSEOPERATIONS_TCC
#define FASTTREE_SSEOPERATIONS_TCC

#include "SSEOperations.h"
#include <emmintrin.h>

#define AbsSSEOperations(...) \
template<typename Precision> \
__VA_ARGS__ fasttree::SSEOperations<Precision>

template<>
void fasttree::SSEOperations<float>::vector_add_mult(float fTot[], float fAdd[], float weight, int64_t n) {
    __m128 w = _mm_set1_ps(weight);
    for (int64_t i = 0; i < n; i += 4) {
        __m128 tot, add;
        tot = _mm_load_ps(fTot + i);
        add = _mm_load_ps(fAdd + i);
        _mm_store_ps(fTot + i, _mm_add_ps(tot, _mm_mul_ps(add, w)));
    }
}

template<>
void fasttree::SSEOperations<double>::vector_add_mult(double fTot[], double fAdd[], double weight, int64_t n) {
    __m128d w = _mm_set1_pd(weight);
    for (int64_t i = 0; i < n; i += 2) {
        __m128d tot, add;
        tot = _mm_load_pd(fTot + i);
        add = _mm_load_pd(fAdd + i);
        _mm_store_pd(fTot + i, _mm_add_pd(tot, _mm_mul_pd(add, w)));
    }
}

#endif

