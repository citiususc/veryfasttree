
#ifndef FASTTREE_SSEOPERATIONS_H
#define FASTTREE_SSEOPERATIONS_H

#include "../DistanceMatrix.h"

namespace fasttree {
    template<typename Precision>
    class SSEOperations {
    public:
        static constexpr int ALIGNMENT = 16;
        typedef Precision numeric_t;

        void vector_add_mult(numeric_t fTot[], numeric_t fAdd[], numeric_t weight, int64_t n);

    };
}

#include "SSEOperations.tcc"

#endif
