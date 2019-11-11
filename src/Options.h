
#ifndef VERYFASTTREE_OPTIONS_H
#define VERYFASTTREE_OPTIONS_H

#include "Constants.h"
#include "Debug.h"
#include <string>
#include <ctime>

namespace veryfasttree {
    struct Options {
        int verbose = 1;
        bool showProgress = true;
        bool slow = false;
        bool fastest = false;
        /* use the second-level top hits heuristic? */
        bool useTopHits2nd = false;
        bool bionj = 0;
        /* 0 means compare nodes to all other nodes */
        double tophitsMult = 1.0;
        /* Parameter for how close is close; also used as a coverage req. */
        double tophitsClose = -1.0;
        /* nTopVisible = m * topvisibleMult; 1 or 2 did not make much difference in either running time or accuracy so
         * I chose a compromise. */
        double topvisibleMult = 1.5;

        /* Refresh if fraction of top-hit-length drops to this */
        double tophitsRefresh = 0.8;
        /* Second-level top heuristic -- only with -fastest */
        double tophits2Mult = 1.0;
        /* Safety factor for second level of top-hits heuristic */
        int tophits2Safety = 3;
        /* Refresh 2nd-level top hits if drops down to this fraction of length */
        double tophits2Refresh = 0.6;
        /* nActive changes by at most this amount before we recompute an out-distance.
         * (Only applies if using the top-hits heuristic) */
        double staleOutLimit = 0.01;
        /* Recompute out profile from scratch if nActive has changed by more than this proportion, and */
        double fResetOutProfile = 0.02;
        /* nActive has also changed more than this amount */
        int nResetOutProfile = 200;
        /* 20 if protein, 4 if nucleotide */
        int nCodes = 20;
        /* If false, use %different as the uncorrected distance */
        bool useMatrix = true;
        /* If true, do a log-correction (scoredist-like or Jukes-Cantor) but only during NNIs and support values,
         * not during neighbor-joining */
        bool logdist = true;
        /* The weight of pseudocounts to avoid artificial long branches when nearby sequences in the tree have little
         * or no overlap (off by default). The prior distance is based on all overlapping positions among the quartet
         * or triplet under consideration. The log correction takes place after the pseudocount is used. */
        double pseudoWeight = 0.0;
        /* Cost of violation of a topological constraint in evolutionary distance or likelihood */
        double constraintWeight = 100.0;
        /* Changes of less than this in tree-length are discounted for purposes of identifying fixed subtrees */
        double MEMinDelta = 1.0e-4;
        bool fastNNI = true;
        /* compute gamma likelihood without reoptimizing branch lengths? */
        bool gammaLogLk = false;
        /* Rounds of optimization of branch lengths; 1 means do 2nd round only if close */
        int mlAccuracy = 1;
        /* Exact or approximate posterior distributions for a.a.s */
        bool exactML = true;

        std::string codesString;/* Protein character */

        int64_t nAlign = 1; /* number of alignments to read */
        std::string matrixPrefix;
        bool makeMatrix = false;
        std::string constraintsFile;
        std::string intreeFile;
        bool intree1 = false;        /* the same starting tree each round */
        int nni = -1;            /* number of rounds of NNI, defaults to 4*log2(n) */
        int spr = 2;            /* number of rounds of SPR */
        int maxSPRLength = 10;    /* maximum distance to move a node */
        int MLnni = -1;        /* number of rounds of ML NNI, defaults to 2*log2(n) */
        bool MLlen = false;        /* optimize branch lengths; no topology changes */
        int nBootstrap = 1000;        /* If set, number of replicates of local bootstrap to do */
        int nRateCats = Constants::nDefaultRateCats;
        bool bUseGtr = false;
        bool bUseLg = false;
        bool bUseWag = false;
        bool bUseGtrRates = false;
        double gtrrates[6];
        bool bUseGtrFreq = false;
        double gtrfreq[4];
        bool bQuote = false;

        std::string inFileName;
        std::string outFileName;
        std::string logFileName;
        bool expert = false;
        long seed = static_cast<long>(std::time(nullptr));

        /*Optimizations*/
        int threads = 1;
        int threadsLevel = 1;
        int threadSubtrees = (unsigned int)~0 >> 1;
        bool doublePrecision = false;
        int fastexp = 0;
        std::string extension;

        /*Precision*/
        double MLMinBranchLengthTolerance;
        double MLFTolBranchLength;
        double MLMinBranchLength;
        double MLMinRelBranchLength;
        double fPostTotalTolerance;

        /*Debug*/
        Debug debug;
    };
}

#endif
