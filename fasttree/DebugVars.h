
#ifndef FASTTREE_DEBUGVARS_H
#define FASTTREE_DEBUGVARS_H

namespace fasttree {

    struct DebugVars {

/* Performance and memory usage */
        long profileOps = 0;        /* Full profile-based distance operations */
        long outprofileOps = 0;        /* How many of profileOps are comparisons to outprofile */
        long seqOps = 0;        /* Faster leaf-based distance operations */
        long profileAvgOps = 0;        /* Number of profile-average steps */
        long nHillBetter = 0;        /* Number of hill-climbing steps */
        long nCloseUsed = 0;        /* Number of "close" neighbors we avoid full search for */
        long nClose2Used = 0;        /* Number of "close" neighbors we use 2nd-level top hits for */
        long nRefreshTopHits = 0;    /* Number of full-blown searches (interior nodes) */
        long nVisibleUpdate = 0;        /* Number of updates of the visible set */
        long nNNI = 0;            /* Number of NNI changes performed */
        long nSPR = 0;            /* Number of SPR changes performed */
        long nML_NNI = 0;        /* Number of max-lik. NNI changes performed */
        long nSuboptimalSplits = 0;    /* # of splits that are rejected given final tree (during bootstrap) */
        long nSuboptimalConstrained = 0; /* Bad splits that are due to constraints */
        long nConstraintViolations = 0;    /* Number of constraint violations */
        long nProfileFreqAlloc = 0;
        long nProfileFreqAvoid = 0;
        long szAllAlloc = 0;
        long mymallocUsed = 0;        /* useful allocations by mymalloc */
        long maxmallocHeap = 0;        /* Maximum of mi.arena+mi.hblkhd from mallinfo (actual mem usage) */
        long nLkCompute = 0;        /* # of likelihood computations for pairs of probability vectors */
        long nPosteriorCompute = 0;    /* # of computations of posterior probabilities */
        long nAAPosteriorExact = 0;    /* # of times compute exact AA posterior */
        long nAAPosteriorRough = 0;    /* # of times use rough approximation */
        long nStarTests = 0;        /* # of times we use star test to avoid testing an NNI */

    };
}

#endif
