
#ifndef FASTTREE_BASICOPERATIONS_TCC
#define FASTTREE_BASICOPERATIONS_TCC

#include "BasicOperations.h"
#include <cmath>

#define AbsBasicOperations(...) \
template<typename Precision> \
__VA_ARGS__ fasttree::BasicOperations<Precision>

AbsBasicOperations(inline void)::vector_multiply(numeric_t f1[], numeric_t f2[], int64_t n, numeric_t fOut[]) {
    for (int64_t i = 0; i < n; i++) {
        fOut[i] = f1[i] * f2[i];
    }
}

AbsBasicOperations(inline Precision)::vector_multiply_sum(numeric_t f1[], numeric_t f2[], int64_t n) {
    numeric_t out = 0.0;
    for (int64_t i = 0; i < n; i++) {
        out += f1[i] * f2[i];
    }
    return out;
}

AbsBasicOperations(inline Precision)::vector_multiply3_sum(numeric_t f1[], numeric_t f2[], numeric_t f3[], int64_t n) {
    numeric_t sum = 0.0;
    for (int64_t i = 0; i < n; i++)
        sum += f1[i] * f2[i] * f3[i];
    return sum;
}

AbsBasicOperations(inline Precision)::vector_dot_product_rot(numeric_t f1[], numeric_t f2[], numeric_t fBy[], int64_t n) {
    numeric_t out1 = 0.0;
    numeric_t out2 = 0.0;
    for (int64_t i = 0; i < n; i++) {
        out1 += f1[i] * fBy[i];
        out2 += f2[i] * fBy[i];
    }
    return out1 * out2;
}

AbsBasicOperations(inline Precision)::vector_sum(numeric_t f1[], int64_t n) {
    numeric_t out = 0.0;
    for (int64_t i = 0; i < n; i++) {
        out += f1[i];
    }
    return (out);
}

AbsBasicOperations(inline void)::vector_multiply_by(numeric_t f[], numeric_t fBy, int64_t n, numeric_t fOut[]) {
    for (int64_t i = 0; i < n; i++) {
        fOut[i] = f[i] * fBy;
    }
}

AbsBasicOperations(inline void)::vector_add_mult(numeric_t fTot[], numeric_t fAdd[], numeric_t weight, int64_t n) {
    for (int64_t i = 0; i < n; i++) {
        fTot[i] += fAdd[i] * weight;
    }
}

AbsBasicOperations(inline void)::matrixt_by_vector4(numeric_t mat[4][4], numeric_t vec[4], numeric_t out[4]) {
    for (int64_t j = 0; j < 4; j++) {
        double sum = 0;
        for (int64_t k = 0; k < 4; k++) {
            sum += vec[k] * mat[k][j];
        }
        out[j] = sum;
    }
}

AbsBasicOperations(inline void)::fastexp(numeric_t fTot[], int64_t n, int lvl) {
    if (lvl == 0) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = (numeric_t) std::exp((double) fTot[k]);
        }
    } else if (lvl == 1) {
        for (int64_t k = 0; k < n; k++) {
            fTot[k] = (numeric_t) std::exp((float) fTot[k]);
        }
    } else if (lvl == 2) {
        int64_t m;
        double xx, px, qx, x;
        _Double u;
        for (int64_t k = 0; k < n; k++) {
            x = (double) fTot[k];
            px = std::floor(1.4426950408889634073599 * x + 0.5);

            m = (int64_t) px;

            x -= px * 6.93145751953125E-1;
            x -= px * 1.42860682030941723212E-6;
            xx = x * x;

            /* px = x * P(x**2). */
            px = 1.26177193074810590878E-4;
            px *= xx;
            px += 3.02994407707441961300E-2;
            px *= xx;
            px += 9.99999999999999999910E-1;
            px *= x;

            /* Evaluate Q(x**2). */
            qx = 3.00198505138664455042E-6;
            qx *= xx;
            qx += 2.52448340349684104192E-3;
            qx *= xx;
            qx += 2.27265548208155028766E-1;
            qx *= xx;
            qx += 2.00000000000000000009E0;

            /* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) ) */
            x = px / (qx - px);
            x = 1.0 + 2.0 * x;

            /* Build 2^n in double. */
            u.d = 0;
            m += 1023;
            u.i = (int64_t) (m) << 52;

            fTot[k] = (numeric_t) x * u.d;
        }
    } else {
        int32_t m;
        float xx, px, qx, x;
        _Float u;
        for (int64_t k = 0; k < n; k++) {
            x = (float) fTot[k];
            px = std::floor(1.4426950408889634073599f * x + 0.5f);

            m = (int32_t) px;

            x -= px * 6.93145751953125E-1f;
            x -= px * 1.42860682030941723212E-6f;
            xx = x * x;

            /* px = x * P(x**2). */
            px = 1.26177193074810590878E-4f;
            px *= xx;
            px += 3.02994407707441961300E-2f;
            px *= xx;
            px += 9.99999999999999999910E-1f;
            px *= x;

            /* Evaluate Q(x**2). */
            qx = 3.00198505138664455042E-6f;
            qx *= xx;
            qx += 2.52448340349684104192E-3f;
            qx *= xx;
            qx += 2.27265548208155028766E-1f;
            qx *= xx;
            qx += 2.00000000000000000009E0f;

            /* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) ) */
            x = px / (qx - px);
            x = 1.0 + 2.0 * x;

            /* Build 2^n in double. */
            u.f = 0;
            m += 127;
            u.i = (int32_t) (m) << 23;

            fTot[k] = (numeric_t) x * u.f;
        }
    }
}

#endif