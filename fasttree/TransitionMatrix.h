
#ifndef FASTTREE_TRANSITIONMATRIX_H
#define FASTTREE_TRANSITIONMATRIX_H

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
namespace fasttree {

    template<typename Precision>
    struct TransitionMatrix {
        typedef Precision numeric_t;

        numeric_t stat[MAXCODES]; /* The stationary distribution */
        numeric_t statinv[MAXCODES];    /* 1/stat */
        /* the eigenmatrix, with the eigenvectors as columns and rotations of individual
           characters as rows. Also includes a NOCODE entry for gaps */
        numeric_t codeFreq[NOCODE + 1][MAXCODES];
        numeric_t eigeninv[MAXCODES][MAXCODES]; /* Inverse of eigenmatrix */
        numeric_t eigeninvT[MAXCODES][MAXCODES]; /* transpose of eigeninv */
        numeric_t eigenval[MAXCODES];    /* Eigenvalues  */
        /* These are for approximate posteriors (off by default) */
        numeric_t nearP[MAXCODES][MAXCODES]; /* nearP[i][j] = P(parent=j | both children are i, both lengths are 0.1 */
        numeric_t nearFreq[MAXCODES][MAXCODES]; /* rotation of nearP/stat */
    };
}

#endif
