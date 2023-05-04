#include "CudaOperations.h"
#include <cublas_v2.h>
#include <omp.h>
#include <exception>
#include <thrust/device_vector.h>
#include <thrust/transform.h>

inline void set_device(int nDevices){
    if(omp_in_parallel()){
        cudaSetDevice(omp_get_thread_num() % nDevices);
    }
}

template<typename Precision>
veryfasttree::CudaOperations<Precision>::CudaOperations(){
    cudaGetDeviceCount(&this->nDevices);
}

template<typename Precision>
void
veryfasttree::CudaOperations<Precision>::vector_multiply(numeric_t f1[], numeric_t f2[], int64_t n, numeric_t fOut[]) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d_1(f1, f1 + n), d_2(f2, f2 + n), d_o(n);
    thrust::transform(d_1.begin(), d_1.end(), d_2.begin(), d_o.begin(), thrust::multiplies<numeric_t>());
    thrust::copy(d_o.begin(), d_o.end(), fOut);

}

template<typename Precision>
Precision veryfasttree::CudaOperations<Precision>::vector_multiply_sum(numeric_t f1[], numeric_t f2[], int64_t n) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d_1(f1, f1 + n), d_2(f2, f2 + n), d_aux(n);
    thrust::transform(d_1.begin(), d_1.end(), d_2.begin(), d_aux.begin(), thrust::multiplies<numeric_t>());
    return thrust::reduce(d_aux.begin(), d_aux.end(), 0, thrust::plus<numeric_t>());
}

template<typename Precision>
Precision veryfasttree::CudaOperations<Precision>::vector_multiply3_sum(numeric_t f1[], numeric_t f2[], numeric_t f3[],
                                                                        int64_t n) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d_1(f1, f1 + n), d_2(f2, f2 + n), d_3(f3, f3 + n), d_aux(n);
    thrust::transform(d_1.begin(), d_1.end(), d_2.begin(), d_aux.begin(), thrust::multiplies<numeric_t>());
    thrust::transform(d_aux.begin(), d_aux.end(), d_3.begin(), d_1.begin(), thrust::multiplies<numeric_t>());
    return thrust::reduce(d_1.begin(), d_1.end(), 0, thrust::plus<numeric_t>());
}


template<typename Precision>
Precision
veryfasttree::CudaOperations<Precision>::vector_dot_product_rot(numeric_t f1[], numeric_t f2[], numeric_t fBy[],
                                                                int64_t n) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d_1(f1, f1 + n), d_2(f2, f2 + n), d_By(fBy, fBy + n), d_out1(n), d_out2(n);
    thrust::transform(d_1.begin(), d_1.end(), d_By.begin(), d_out1.begin(), thrust::multiplies<numeric_t>());
    thrust::transform(d_2.begin(), d_2.end(), d_By.begin(), d_out2.begin(), thrust::multiplies<numeric_t>());
    return thrust::reduce(d_out1.begin(), d_out1.end(), 0, thrust::plus<numeric_t>()) *
           thrust::reduce(d_out2.begin(), d_out2.end(), 0, thrust::plus<numeric_t>());
}


template<typename Precision>
void veryfasttree::CudaOperations<Precision>::vector_add(numeric_t fTot[], numeric_t fAdd[], int64_t n) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d_Tot(fTot, fTot + n), d_Add(fAdd, fAdd + n), d_o(n);
    thrust::transform(d_Tot.begin(), d_Tot.end(), d_Add.begin(), d_o.begin(), thrust::plus<numeric_t>());
    thrust::copy(d_o.begin(), d_o.end(), fTot);
}

template<typename Precision>
Precision veryfasttree::CudaOperations<Precision>::vector_sum(numeric_t f1[], int64_t n) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d_1(f1, f1 + n);
    return thrust::reduce(d_1.begin(), d_1.end(), 0, thrust::plus<numeric_t>());
}

template<typename Precision>
void veryfasttree::CudaOperations<Precision>::vector_multiply_by(numeric_t f[], numeric_t fBy, int64_t n,
                                                                 numeric_t fOut[]) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d(f, f + n), d_By(n, fBy), d_o(n);
    thrust::transform(d.begin(), d.end(), d_By.begin(), d_o.begin(), thrust::multiplies<numeric_t>());
    thrust::copy(d_o.begin(), d_o.end(), fOut);
}

template<typename Precision>
void veryfasttree::CudaOperations<Precision>::vector_add_mult(numeric_t fTot[], numeric_t fAdd[], numeric_t weight,
                                                              int64_t n) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d_Tot(fTot, fTot + n), d_Add(fAdd, fAdd + n), d_weight(n, weight), d_o(n);
    thrust::transform(d_Add.begin(), d_Add.end(), d_weight.begin(), d_o.begin(), thrust::multiplies<numeric_t>());
    thrust::transform(d_Tot.begin(), d_Tot.end(), d_o.begin(), d_Add.begin(), thrust::plus<numeric_t>());
    thrust::copy(d_Add.begin(), d_Add.end(), fTot);
}

template<>
template<int row>
void veryfasttree::CudaOperations<float>::
matrix_by_vector4(numeric_t mat[][row], numeric_t vec[], numeric_t out[]) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d_A((numeric_t*)mat, (numeric_t*)mat + 16), d_vec(vec, vec + 4), d_o(4);
    const float alpha = 1;
    const float beta = 0;
    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 4, 1, 4,
                &alpha, thrust::raw_pointer_cast(&d_A[0]), 4, thrust::raw_pointer_cast(&d_vec[0]), 4,
                &beta, thrust::raw_pointer_cast(&d_o[0]), 4);
    cublasDestroy(handle);
    thrust::copy(d_o.begin(), d_o.end(), out);
}

template<>
template<int row>
void veryfasttree::CudaOperations<double>::
matrix_by_vector4(numeric_t mat[][row], numeric_t vec[], numeric_t out[]) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d_A((numeric_t*)mat, (numeric_t*)mat + 16), d_vec(vec, vec + 4), d_o(4);
    const double alpha = 1;
    const double beta = 0;
    cublasHandle_t handle;
    cublasCreate(&handle);
    cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, 4, 1, 4,
                &alpha, thrust::raw_pointer_cast(&d_A[0]), 4, thrust::raw_pointer_cast(&d_vec[0]), 4,
                &beta, thrust::raw_pointer_cast(&d_o[0]), 4);
    cublasDestroy(handle);
    thrust::copy(d_o.begin(), d_o.end(), out);
}

struct Expf {
    __device__ float operator()(const float &x) const {
        return expf(x);
    }
};

template<>
void veryfasttree::CudaOperations<float>::fastexp(numeric_t fTot[], int64_t n, int lvl) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d_Tot(fTot, fTot + n), d_o(n);
    thrust::transform(d_Tot.begin(), d_Tot.end(), d_o.begin(), Expf());
    thrust::copy(d_o.begin(), d_o.end(), fTot);
}

struct Expd {
    __device__ double operator()(const double &x) const {
        return exp(x);
    }
};

template<>
void veryfasttree::CudaOperations<double>::fastexp(numeric_t fTot[], int64_t n, int lvl) {
    set_device(this->nDevices);
    thrust::device_vector <numeric_t> d_Tot(fTot, fTot + n), d_o(n);
    thrust::transform(d_Tot.begin(), d_Tot.end(), d_o.begin(), Expd());
    thrust::copy(d_o.begin(), d_o.end(), fTot);
}

template
class veryfasttree::CudaOperations<float>;

template void veryfasttree::CudaOperations<float>::matrix_by_vector4(float mat[][4], float vec[], float out[]);
template void veryfasttree::CudaOperations<float>::matrix_by_vector4(float mat[][20], float vec[], float out[]);

template
class veryfasttree::CudaOperations<double>;

template void veryfasttree::CudaOperations<double>::matrix_by_vector4(double mat[][4], double vec[], double out[]);
template void veryfasttree::CudaOperations<double>::matrix_by_vector4(double mat[][20], double vec[], double out[]);

void veryfasttree::configCuda(Options& options){
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    if(nDevices == 0){
        throw std::runtime_error("No CUDA compatible device found");
    }
}


