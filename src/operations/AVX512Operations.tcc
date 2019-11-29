
#ifndef FASTTREE_AVX512OPERATIONS_TCC
#define FASTTREE_AVX512OPERATIONS_TCC

#include "AVX512Operations.tcc"

template<>
inline void veryfasttree::AVX512Operations<float>::vector_multiply(float f1[], float f2[], int64_t n, float fOut[]) {

}

template<>
inline void veryfasttree::AVX512Operations<double>::vector_multiply(double f1[], double f2[], int64_t n, double fOut[]) {

}

template<>
inline float veryfasttree::AVX512Operations<float>::vector_multiply_sum(float f1[], float f2[], int64_t n) {
    return 0;
}

template<>
inline double veryfasttree::AVX512Operations<double>::vector_multiply_sum(double f1[], double f2[], int64_t n) {
    return 0;
}

template<>
inline float
veryfasttree::AVX512Operations<float>::vector_multiply3_sum(float f1[], float f2[], float f3[], int64_t n) {
    return 0;
}

template<>
inline double
veryfasttree::AVX512Operations<double>::vector_multiply3_sum(double f1[], double f2[], double f3[], int64_t n) {
    return 0;
}

template<>
inline float
veryfasttree::AVX512Operations<float>::vector_dot_product_rot(float f1[], float f2[], float fBy[], int64_t n) {
    return 0;
}

template<>
inline double veryfasttree::AVX512Operations<double>::
vector_dot_product_rot(double f1[], double f2[], double fBy[], int64_t n) {
    return 0;
}

inline void veryfasttree::AVX512Operations<float>::vector_add(float fTot[], float fAdd[], int64_t n){

}

inline void veryfasttree::AVX512Operations<double>::vector_add(double fTot[], double fAdd[], int64_t n){

}

template<>
inline float veryfasttree::AVX512Operations<float>::vector_sum(float f1[], int64_t n) {
    return 0;
}

template<>
inline double veryfasttree::AVX512Operations<double>::vector_sum(double f1[], int64_t n) {
    return 0;
}

template<>
inline void veryfasttree::AVX512Operations<float>::vector_multiply_by(float f[], float fBy, int64_t n, float fOut[]) {

}

template<>
inline void veryfasttree::AVX512Operations<double>::vector_multiply_by(double f[], double fBy, int64_t n, double fOut[]) {

}

template<>
inline void veryfasttree::AVX512Operations<float>::vector_add_mult(float fTot[], float fAdd[], float weight, int64_t n) {

}

template<>
inline void
veryfasttree::AVX512Operations<double>::vector_add_mult(double fTot[], double fAdd[], double weight, int64_t n) {

}

template<>
template<int row>
inline void veryfasttree::AVX512Operations<float>::matrixt_by_vector4(float mat[][row], float vec[], float out[]) {

}

template<>
template<int row>
inline void veryfasttree::AVX512Operations<double>::matrixt_by_vector4(double mat[][row], double vec[], double out[]) {

}

template<>
inline void veryfasttree::AVX512Operations<float>::fastexp(float fTot[], int64_t n, int lvl) {

}

template<>
inline void veryfasttree::AVX512Operations<double>::fastexp(double fTot[], int64_t n, int lvl) {

}

#endif
