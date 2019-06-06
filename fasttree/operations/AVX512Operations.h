
#ifndef FASTTREE_AVX512OPERATIONS_H
#define FASTTREE_AVX512OPERATIONS_H

#include <boost/align/aligned_allocator.hpp>
#include <immintrin.h>

namespace fasttree {
    template<typename Precision>
    class AVX512Operations {
    public:
        static constexpr int ALIGNMENT = 64;
        using Allocator = boost::alignment::aligned_allocator<Precision, ALIGNMENT>;
        typedef Precision numeric_t;

        void vector_multiply(numeric_t f1[], numeric_t f2[], int64_t n, numeric_t fOut[]);

        numeric_t vector_multiply_sum(numeric_t f1[], numeric_t f2[], int64_t n);

        numeric_t vector_multiply3_sum(numeric_t f1[], numeric_t f2[], numeric_t f3[], int64_t n);

        numeric_t vector_dot_product_rot(numeric_t f1[], numeric_t f2[], numeric_t fBy[], int64_t n);

        numeric_t vector_sum(numeric_t f1[], int64_t n);

        void vector_multiply_by(numeric_t f[], numeric_t fBy, int64_t n, numeric_t fOut[]);

        void vector_add_mult(numeric_t fTot[], numeric_t fAdd[], numeric_t weight, int64_t n);

        void matrixt_by_vector4(numeric_t mat[4][4], numeric_t vec[4], numeric_t out[4]);

        void fastexp(numeric_t fTot[], int64_t n, int lvl);

    };
}

#include "AVX512Operations.tcc"

#endif
