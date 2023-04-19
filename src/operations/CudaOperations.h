

#ifndef VERYFASTTREE_CUDAOPERATION_H
#define VERYFASTTREE_CUDAOPERATION_H

#include <memory>
#include "../Options.h"

namespace veryfasttree {

    template<typename Precision>
    class CudaOperations {
    public:        /*
         * No alignment required, values less than sizeof(Precision) are ignored,
         * value 0 is not supported by gcc 5 and 6.
         */
        static constexpr int ALIGNMENT = sizeof(Precision);
        using Allocator = std::allocator<Precision>;
        typedef Precision numeric_t;
        void vector_multiply(numeric_t f1[], numeric_t f2[], int64_t n, numeric_t fOut[]);

        numeric_t vector_multiply_sum(numeric_t f1[], numeric_t f2[], int64_t n);

        numeric_t vector_multiply3_sum(numeric_t f1[], numeric_t f2[], numeric_t f3[], int64_t n);

        numeric_t vector_dot_product_rot(numeric_t f1[], numeric_t f2[], numeric_t fBy[], int64_t n);

        void vector_add(numeric_t fTot[], numeric_t fAdd[], int64_t n);

        numeric_t vector_sum(numeric_t f1[], int64_t n);

        void vector_multiply_by(numeric_t f[], numeric_t fBy, int64_t n, numeric_t fOut[]);

        void vector_add_mult(numeric_t fTot[], numeric_t fAdd[], numeric_t weight, int64_t n);

        template <int row>
        void matrix_by_vector4(numeric_t mat[][row], numeric_t vec[], numeric_t out[]);

        void fastexp(numeric_t fTot[], int64_t n, int lvl);

        CudaOperations();
    private:
        int nDevices;
    };

    void configCuda(Options& options);
}

#endif //VERYFASTTREE_CUDAOPERATION_H
