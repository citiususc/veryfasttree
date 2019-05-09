
#ifndef FASTTREE_BASICOPERATIONS_TCC
#define FASTTREE_BASICOPERATIONS_TCC

#include "BasicOperations.h"

#define AbsBasicOperations(...) \
template<typename Precision> \
__VA_ARGS__ fasttree::BasicOperations<Precision>

AbsBasicOperations(void)::vector_multiply(numeric_t f1[], numeric_t f2[], int64_t n, numeric_t fOut[]){
    for (int64_t i = 0; i < n; i++) {
        fOut[i] = f1[i] * f2[i];
    }
}

AbsBasicOperations(Precision)::vector_multiply_sum(numeric_t f1[], numeric_t f2[], int64_t n){
    numeric_t out = 0.0;
    for (int64_t i=0; i < n; i++) {
        out += f1[i] * f2[i];
    }
    return out;
}

AbsBasicOperations(Precision)::vector_multiply3_sum(numeric_t f1[], numeric_t f2[], numeric_t f3[], int64_t n){
    numeric_t sum = 0.0;
    for (int64_t i = 0; i < n; i++)
        sum += f1[i]*f2[i]*f3[i];
    return sum;
}

AbsBasicOperations(Precision)::vector_dot_product_rot(numeric_t f1[], numeric_t f2[], numeric_t fBy[], int64_t n){
    numeric_t out1 = 0.0;
    numeric_t out2 = 0.0;
    for (int64_t i=0; i < n; i++) {
        out1 += f1[i]*fBy[i];
        out2 += f2[i]*fBy[i];
    }
    return out1*out2;
}

AbsBasicOperations(Precision)::vector_sum(numeric_t f1[], int64_t n){
    numeric_t out = 0.0;
    for (int64_t i = 0; i < n; i++) {
        out += f1[i];
    }
    return(out);
}

AbsBasicOperations(void)::vector_multiply_by(numeric_t f[], numeric_t fBy, int64_t n){
    for (int64_t i = 0; i < n; i++) {
        f[i] *= fBy;
    }
}

AbsBasicOperations(void)::vector_add_mult(numeric_t fTot[], numeric_t fAdd[], numeric_t weight, int64_t n) {
    for (int64_t i = 0; i < n; i++) {
        fTot[i] += fAdd[i] * weight;
    }
}

AbsBasicOperations(void)::matrixt_by_vector4(numeric_t mat[4][MAXCODES], numeric_t vec[4], numeric_t out[4]){
    for (int64_t j = 0; j < 4; j++) {
        double sum = 0;
        for (int64_t k = 0; k < 4; k++) {
            sum += vec[k] * mat[k][j];
        }
        out[j] = sum;
    }
}

#endif