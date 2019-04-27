
#ifndef FASTTREE_NEIGHBOURJOINING_H
#define FASTTREE_NEIGHBOURJOINING_H

#include "DistanceMatrix.h"
#include "TransitionMatrix.h"
#include <vector>
#include <iostream>


namespace fasttree {

    template<typename Precision>
    class NeighbourJoining {
    public:
        typedef Precision numeric_t;

        struct Profile {
            /* alignment profile */
            std::vector<numeric_t> weights;
            std::string codes;
            std::vector<numeric_t> vectors;        /* empty if no non-constant positions, e.g. for leaves */
            std::vector<numeric_t> codeDist;        /* Optional -- distance to each code at each position */

            /* constraint profile */
            std::vector<bool> nOn;
            std::vector<bool> nOff;

            numeric_t nGaps; /*precalculated in construction*/

            Profile(size_t nPos, size_t nConstraints);

            Profile();
        };

        struct Rates {
            std::vector<numeric_t> rates;    /* 1 per rate category */
            std::vector<size_t>  ratecat;    /* 1 category per position */

            /* Allocate or reallocate the rate categories, and set every position
               to category 0 and every category's rate to 1.0
               If nRateCategories=0, just deallocate
            */
            Rates(size_t nRateCategories, size_t nPos);
        };

        struct Children {
            int nChild = 0;
            size_t child[3];
        };

        NeighbourJoining(const Options &options, std::ostream &log, const std::vector<std::string> &seqs, size_t nPos,
                         const std::vector<std::string> &constraintSeqs,
                         DistanceMatrix<Precision> &distanceMatrix, TransitionMatrix<Precision> &transmat);

    private:
        std::ostream &log;
        const Options &options;
        /* The input */
        size_t nPos;
        const std::vector<std::string> &seqs;    /* the aligment sequences array (not reallocated) */
        DistanceMatrix<Precision> &distanceMatrix; /* a ref, not set if using %identity distance */
        TransitionMatrix<Precision> &transmat; /* a ref, not set for Jukes-Cantor */
        /* Topological constraints are represented for each sequence as binary characters
           with values of '0', '1', or '-' (for missing data)
           Sequences that have no constraint may have a NULL string
        */
        const std::vector<std::string> &constraintSeqs;

        /* The profile data structures */
        size_t maxnode;            /* The next index to allocate */
        size_t maxnodes;            /* Space allocated in data structures below */
        std::vector<Profile> profiles; /* Profiles of leaves and intermediate nodes */
        std::vector<numeric_t> diameter;        /* To correct for distance "up" from children (if any) */
        std::vector<numeric_t> varDiameter;        /* To correct variances for distance "up" */
        std::vector<numeric_t> selfdist;        /* Saved for use in some formulas */
        std::vector<numeric_t> selfweight;        /* Saved for use in some formulas */

        /* Average profile of all active nodes, the "outprofile"
         * If all inputs are ungapped, this has weight 1 (not nSequences) at each position
         * The frequencies all sum to one (or that is implied by the eigen-representation)
         */
        Profile outprofile;
        double totdiam;

        /* We sometimes use stale out-distances, so we remember what nActive was  */
        std::vector<numeric_t> outDistances;        /* Sum of distances to other active (parent==-1) nodes */
        std::vector<size_t> nOutDistActive;        /* What nActive was when this outDistance was computed */

        /* the inferred tree */
        int64_t root;                   /* index of the root. Unlike other internal nodes, it has 3 children */
        std::vector<int64_t> parent;    /* -1 or index of parent */
        std::vector<Children> child;
        std::vector<numeric_t> branchlength;        /* Distance to parent */
        std::vector<numeric_t> support;        /* 1 for high-confidence nodes */

        /* auxilliary data for maximum likelihood (defaults to 1 category of rate=1.0) */
        Rates rates;

        uint8_t charToCode[256];

        void
        seqToProfile(size_t i, Profile &profile, const std::string &seq, const std::string &constraintSeq);

        /* outProfile() computes the out-profile */
        void outProfile(Profile &out);

        /* Use out-profile and NJ->totdiam to recompute out-distance for node iNode
           Only does this computation if the out-distance is "stale" (nOutDistActive[iNode] != nActive)
           Note "IN/UPDATE" for NJ always means that we may update out-distances but otherwise
           make no changes.
        */
        void setOutDistance(size_t i, size_t nActive);

        /* only handles leaf sequences */
        size_t nGaps(size_t i);
    };
}

#include "NeighbourJoining.tcc"

#endif
