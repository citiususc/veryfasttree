
#ifndef FASTTREE_SSE128OPERATIONS_H
#define FASTTREE_SSE128OPERATIONS_H

#include <boost/align/aligned_allocator.hpp>
#include <immintrin.h>

namespace fasttree {
    template<typename Precision>
    class SSE128Operations {
    public:
        static constexpr int ALIGNMENT = 16;
        using Allocator = boost::alignment::aligned_allocator<Precision, ALIGNMENT>;
        typedef Precision numeric_t;

        inline void vector_multiply(numeric_t f1[], numeric_t f2[], int64_t n, numeric_t fOut[]);

        inline numeric_t vector_multiply_sum(numeric_t f1[], numeric_t f2[], int64_t n);

        inline numeric_t vector_multiply3_sum(numeric_t f1[], numeric_t f2[], numeric_t f3[], int64_t n);

        inline numeric_t vector_dot_product_rot(numeric_t f1[], numeric_t f2[], numeric_t fBy[], int64_t n);

        inline numeric_t vector_sum(numeric_t f1[], int64_t n);

        inline void vector_multiply_by(numeric_t f[], numeric_t fBy, int64_t n, numeric_t fOut[]);

        inline void vector_add_mult(numeric_t fTot[], numeric_t fAdd[], numeric_t weight, int64_t n);

        inline void matrixt_by_vector4(numeric_t mat[4][4], numeric_t vec[4], numeric_t out[4]);

        inline void fastexp(numeric_t fTot[], int64_t n, int lvl);

    private:
        template<typename Tp>
        inline numeric_t mm_sum(Tp sum);

        inline __m128 fastexpImpl(__m128 vx);

        inline __m128d fastexpImpl(__m128d vx);
    };
}

/*
 * A template specialization must be declared inside namespace in gcc 5 and 6.
 */
namespace fasttree {

#include "SSE128Operations.tcc"

}

#endif
