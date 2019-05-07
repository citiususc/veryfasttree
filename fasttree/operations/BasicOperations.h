
#ifndef FASTTREE_BASICOPERATIONS_H
#define FASTTREE_BASICOPERATIONS_H

#include "../DistanceMatrix.h"

namespace fasttree {

    template<typename Precision>
    class BasicOperations {
    public:
        typedef Precision numeric_t;

        numeric_t vector_multiply_sum(numeric_t f1[], numeric_t f2[], int64_t n);

        numeric_t vector_multiply3_sum(numeric_t f1[], numeric_t f2[], numeric_t f3[], int64_t n);

        numeric_t vector_dot_product_rot(numeric_t f1[], numeric_t f2[], numeric_t fBy[], int64_t n);

        numeric_t vector_sum(numeric_t f1[], int64_t n);

        void vector_multiply_by(numeric_t f[], numeric_t fBy, int64_t n);

        void vector_add_mult(numeric_t fTot[], numeric_t fAdd[], numeric_t weight, int64_t n);

        void matrixt_by_vector4(numeric_t mat[4][MAXCODES], numeric_t vec[4], numeric_t out[4]);

    };
}

#include "BasicOperations.tcc"

#endif
