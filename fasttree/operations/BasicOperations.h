
#ifndef FASTTREE_BASICOPERATIONS_H
#define FASTTREE_BASICOPERATIONS_H


namespace fasttree {

    template<typename Precision>
    class BasicOperations {
    public:
        /*
         * No alignment required, values less than sizeof(Precision) are ignored,
         * value 0 is not supported by gcc 5 and 6.
         */
        static constexpr int ALIGNMENT = 2;
        using Allocator = std::allocator<Precision>;
        typedef Precision numeric_t;

        inline void vector_multiply(numeric_t f1[], numeric_t f2[], int64_t n, numeric_t fOut[]);

        inline numeric_t vector_multiply_sum(numeric_t f1[], numeric_t f2[], int64_t n);

        inline numeric_t vector_multiply3_sum(numeric_t f1[], numeric_t f2[], numeric_t f3[], int64_t n);

        inline numeric_t vector_dot_product_rot(numeric_t f1[], numeric_t f2[], numeric_t fBy[], int64_t n);

        inline numeric_t vector_sum(numeric_t f1[], int64_t n);

        inline void vector_multiply_by(numeric_t f[], numeric_t fBy, int64_t n, numeric_t fOut[]);

        inline void vector_add_mult(numeric_t fTot[], numeric_t fAdd[], numeric_t weight, int64_t n);

        template <int row>
        inline void matrixt_by_vector4(numeric_t mat[][row], numeric_t vec[], numeric_t out[]);

        inline void fastexp(numeric_t fTot[], int64_t n, int lvl);

    private:
        union _Float{
            float f;
            int32_t i;
        };

        union _Double{
            double d;
            int64_t i;
        };

    };
}

#include "BasicOperations.tcc"

#endif
