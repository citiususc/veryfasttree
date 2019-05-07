
#ifndef FASTTREE_DEBUG_H
#define FASTTREE_DEBUG_H

#include <memory>

namespace fasttree {

    struct Debug {

        /* Performance and memory usage */
        int64_t profileOps = 0;        /* Full profile-based distance operations */
        int64_t outprofileOps = 0;        /* How many of profileOps are comparisons to outprofile */
        int64_t seqOps = 0;        /* Faster leaf-based distance operations */
        int64_t profileAvgOps = 0;        /* Number of profile-average steps */
        int64_t nHillBetter = 0;        /* Number of hill-climbing steps */
        int64_t nCloseUsed = 0;        /* Number of "close" neighbors we avoid full search for */
        int64_t nClose2Used = 0;        /* Number of "close" neighbors we use 2nd-level top hits for */
        int64_t nRefreshTopHits = 0;    /* Number of full-blown searches (interior nodes) */
        int64_t nVisibleUpdate = 0;        /* Number of updates of the visible set */
        int64_t nNNI = 0;            /* Number of NNI changes performed */
        int64_t nSPR = 0;            /* Number of SPR changes performed */
        int64_t nML_NNI = 0;        /* Number of max-lik. NNI changes performed */
        int64_t nSuboptimalSplits = 0;    /* # of splits that are rejected given final tree (during bootstrap) */
        int64_t nSuboptimalConstrained = 0; /* Bad splits that are due to constraints */
        int64_t nConstraintViolations = 0;    /* Number of constraint violations */
        int64_t nProfileFreqAlloc = 0;
        int64_t nProfileFreqAvoid = 0;
        int64_t szAllAlloc = 0;
        int64_t mymallocUsed = 0;        /* useful allocations by mymalloc */
        int64_t maxmallocHeap = 0;        /* Maximum of mi.arena+mi.hblkhd from mallinfo (actual mem usage) */
        int64_t nLkCompute = 0;        /* # of likelihood computations for pairs of probability vectors */
        int64_t nPosteriorCompute = 0;    /* # of computations of posterior probabilities */
        int64_t nAAPosteriorExact = 0;    /* # of times compute exact AA posterior */
        int64_t nAAPosteriorRough = 0;    /* # of times use rough approximation */
        int64_t nStarTests = 0;        /* # of times we use star test to avoid testing an NNI */

        void reset() {
            (*this) = Debug();
        }
    };

}

#endif
