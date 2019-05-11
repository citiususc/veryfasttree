
#ifndef FASTTREE_CONSTANTS_H
#define FASTTREE_CONSTANTS_H

#include <string>

namespace fasttree {

    namespace Constants {

        /* Maximum likelihood options and constants */
        /* These are used to rescale likelihood values and avoid taking a logarithm at each position */
        const double LkUnderflow = 1.0e-4;
        const double LkUnderflowInv = 1.0e4;
        const double LogLkUnderflow = 9.21034037197618; /* -log(LkUnderflowInv) */
        const double Log2 = 0.693147180559945;
        /* These are used to limit the optimization of branch lengths.
           Also very short branch lengths can create numerical problems.
           In version 2.1.7, the minimum branch lengths (MLMinBranchLength and MLMinRelBranchLength)
           were increased to prevent numerical problems in rare cases.
           In version 2.1.8, to provide useful branch lengths for genome-wide alignments,
           the minimum branch lengths were dramatically decreased if USE_DOUBLE is defined.
        */

        /* absolute tolerance for optimizing branch lengths */
        const double MLMinBranchLengthTolerance = 1.0e-9;
        const double MLMinBranchLengthToleranceDouble = 1.0e-4;
        /* fractional tolerance for optimizing branch lengths */
        const double MLFTolBranchLength = 0.001;
        const double MLFTolBranchLengthDouble = 0.001;
        /* minimum value for branch length */
        const double MLMinBranchLength = 5.0e-9;
        const double MLMinBranchLengthDouble = 5.0e-4;
        /* minimum of rate * length */
        const double MLMinRelBranchLength = 2.5e-9;
        const double MLMinRelBranchLengthDouble = 2.5e-4;
        /* posterior vector must sum to at least this before rescaling */
        const double fPostTotalTolerance = 1.0e-20;
        const double fPostTotalToleranceDouble = 1.0e-10;

        const double closeLogLkLimit = 5.0;    /* If partial optimization of an NNI looks like it would decrease the log likelihood
				   by this much or more then do not optimize it further */
        const double treeLogLkDelta = 0.1;    /* Give up if tree log-lk changes by less than this; NNIs that change
				   likelihood by less than this also are considered unimportant
				   by some heuristics */
        const double approxMLminf = 0.95;    /* Only try to approximate posterior distributions if max. value is at least this high */
        const double approxMLminratio =
                2 / 3.0;/* Ratio of approximated/true posterior values must be at least this high */
        const double approxMLnearT = 0.2;    /* 2nd component of near-constant posterior distribution uses this time scale */
        const int nDefaultRateCats = 20;

        /* Protein character set */
        const std::string codesStringAA = "ARNDCQEGHILKMFPSTWYV";
        const std::string codesStringNT = "ACGT";

        const std::string version = "x.x.x";
        const std::string compileFlags =
                "(OpenMP"
                #ifdef USE_SSE3
                ", SSE3"
                #endif
                #ifdef USE_AVX
                ", AVX"
                #endif
                #ifdef USE_AVX2
                ", AVX2"
                #endif
                #ifdef USE_AVX512
                ", AVX512"
                #endif
                ")";
    };

}

#endif
