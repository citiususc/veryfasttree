
#ifndef VERYFASTTREE_TRANSITIONMATRIX_H
#define VERYFASTTREE_TRANSITIONMATRIX_H

#include "DistanceMatrix.h"

#define NOCODE 127

/* A transition matrix gives the instantaneous rate of change of frequencies
   df/dt = M . f
   which is solved by
   f(t) = exp(M) . f(0)
   and which is not a symmetric matrix because of
   non-uniform stationary frequencies stat, so that
   M stat = 0
   M(i,j) is instantaneous rate of j -> i, not of i -> j

   S = diag(sqrt(stat)) is a correction so that
   M' = S**-1 M S is symmetric
   Let W L W**-1 = M' be an eigendecomposition of M'
   Because M' is symmetric, W can be a rotation, and W**-1 = t(W)
   Set V = S*W
   M = V L V**-1 is an eigendecomposition of M
   Note V**-1 = W**-1 S**-1 = t(W) S**-1

   Evolution by time t is given by

   exp(M*t) = V exp(L*t) V**-1
   P(A & B | t) = B . exp(M*t) . (A * stat)
   note this is *not* the same as P(A->B | t)

   and we can reduce some of the computations from O(a**2) to O(a) time,
   where a is the alphabet size, by storing frequency vectors as
   t(V) . f = t(W) . t(S) . f

   Then
   P(f0 & f1 | t) = f1 . exp(M*t) . f0 * (f0 . stat) = sum(r0j * r1j * exp(l_j*t))
   where r0 and r1 are the transformed vectors

   Posterior distribution of P given children f0 and f1 is given by
   P(i | f0, f1, t0, t1) = stat * P(i->f0 | t0) * P(i->f1 | t1)
   = P(i & f0 | t0) * P(i & f1 | t1) / stat
   ~ (V . exp(t0*L) . r0) * (V . exp(t1*L) . r1) / stat

   When normalize this posterior distribution (to sum to 1), divide by stat,
   and transform by t(V) -- this is the "profile" of internal nodes

   To eliminate the O(N**2) step of transforming by t(V), if the posterior
   distribution of an amino acid is near 1 then we can approximate it by
   P(i) ~= (i==A) * w + nearP(i) * (1-w), where
   w is fit so that P(i==A) is correct
   nearP = Posterior(i | i, i, 0.1, 0.1) [0.1 is an arbitrary choice]
   and we confirm that the approximation works well before we use it.

   Given this parameter w we can set
   rotated_posterior = rotation(w * (i==A)/stat + (1-w) * nearP/stat)
   = codeFreq(A) * w/stat(A) + nearFreq(A) * (1-w)
 */
namespace veryfasttree {

    template<typename Precision, int Aligment>
    struct TransitionMatrix {
        typedef Precision numeric_t;

        alignas(Aligment) numeric_t stat[MAXCODES]; /* The stationary distribution */
        alignas(Aligment) numeric_t statinv[MAXCODES];    /* 1/stat */
        /* the eigenmatrix, with the eigenvectors as columns and rotations of individual
           characters as rows. Also includes a NOCODE entry for gaps */
        alignas(Aligment) numeric_t codeFreq[NOCODE + 1][alignsz(MAXCODES, Aligment / sizeof(numeric_t))];
        /* Inverse of eigenmatrix */
        alignas(Aligment) numeric_t eigeninv[MAXCODES][alignsz(MAXCODES, Aligment / sizeof(numeric_t))];
        alignas(Aligment) numeric_t eigeninvT[MAXCODES][MAXCODES]; /* transpose of eigeninv */
        alignas(Aligment) numeric_t eigenval[MAXCODES];    /* Eigenvalues  */
        /* These are for approximate posteriors (off by default) */
        alignas(Aligment) numeric_t nearP[MAXCODES][MAXCODES]; /* nearP[i][j] = P(parent=j | both children are i, both lengths are 0.1 */
        alignas(Aligment) numeric_t nearFreq[MAXCODES][MAXCODES]; /* rotation of nearP/stat */

        TransitionMatrix();

        void createTransitionMatrixJTT92(const Options &options);

        void createTransitionMatrixWAG01(const Options &options);

        void createTransitionMatrixLG08(const Options &options);

        void createGTR(const Options &options, double gtrrates[]/*ac,ag,at,cg,ct,gt*/, double gtrfreq[]/*ACGT*/);

        void readAATransitionMatrix(const Options &options, /*IN*/ const std::string& filename);

        operator bool();

    private:

        bool setted;

        /* Takes as input the transpose of the matrix V, with i -> j
           This routine takes care of setting the diagonals
        */
        void createTransitionMatrix(const Options &options, const double matrix[MAXCODES][MAXCODES],
                                    const double stat[MAXCODES]);

        /* Numerical recipes code for eigen decomposition (actually taken from RAxML rev_functions.c) */
        void tred2(double a[], const int n, const int np, double d[], double e[]);

        void tqli(double d[], double e[], int n, int np, double z[]);

        inline double pythag(double a, double b);

        static const double statJTT92[MAXCODES];
        static const double matrixJTT92[MAXCODES][MAXCODES];
        static const double statWAG01[MAXCODES];
        static const double matrixWAG01[MAXCODES][MAXCODES];
        static const double statLG08[MAXCODES];
        static const double matrixLG08[MAXCODES][MAXCODES];
    };
}

#include "TransitionMatrix.tcc"

#endif
