
#ifndef VERYFASTTREE_AVX512OPERATIONS_H
#define VERYFASTTREE_AVX512OPERATIONS_H

#include <boost/align/aligned_allocator.hpp>
#include <immintrin.h>

namespace veryfasttree {
    template<typename Precision>
    class AVX512Operations {
    public:
        static constexpr int ALIGNMENT = 64;
        using Allocator = boost::alignment::aligned_allocator<Precision, ALIGNMENT>;
        typedef Precision numeric_t;

        inline void vector_multiply(numeric_t f1[], numeric_t f2[], int64_t n, numeric_t fOut[]);

        inline numeric_t vector_multiply_sum(numeric_t f1[], numeric_t f2[], int64_t n);

        inline numeric_t vector_multiply3_sum(numeric_t f1[], numeric_t f2[], numeric_t f3[], int64_t n);

        inline numeric_t vector_dot_product_rot(numeric_t f1[], numeric_t f2[], numeric_t fBy[], int64_t n);

        inline void vector_add(numeric_t fTot[], numeric_t fAdd[], int64_t n);

        inline numeric_t vector_sum(numeric_t f1[], int64_t n);

        inline void vector_multiply_by(numeric_t f[], numeric_t fBy, int64_t n, numeric_t fOut[]);

        inline void vector_add_mult(numeric_t fTot[], numeric_t fAdd[], numeric_t weight, int64_t n);

        template <int row>
        inline void matrixt_by_vector4(numeric_t mat[][row], numeric_t vec[], numeric_t out[]);

        inline void fastexp(numeric_t fTot[], int64_t n, int lvl);
    private:

        inline numeric_t mm_sum(__m128 sum);

        inline numeric_t mm_sum(__m256d sum);

        inline __m128 fastexpImpl(__m128 vx);

        inline __m256d fastexpImpl(__m256d vx);

        inline __m512 fastexpImpl(__m512 vx);

        inline __m512d fastexpImpl(__m512d vx);

    };
}

/*
 * A template specialization must be declared inside namespace in gcc 5 and 6.
 */
namespace veryfasttree {

#include "AVX512Operations.tcc"

}
#endif
