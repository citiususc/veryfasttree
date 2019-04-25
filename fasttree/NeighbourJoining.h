
#ifndef FASTTREE_NEIGHBOURJOINING_H
#define FASTTREE_NEIGHBOURJOINING_H

namespace fasttree {

    template<typename Precision>
    class NeighbourJoining {
    public:
        typedef Precision numeric_t;

    private:
        /* The input */
        int nSeq;
        int nPos;
        char **seqs;			/* the aligment sequences array (not reallocated) */
        //distance_matrix_t *distance_matrix; /* a pointer (not reallocated), or NULL if using %identity distance */
        //transition_matrix_t *transmat; /* a pointer (is allocated), or NULL for Jukes-Cantor */
        /* Topological constraints are represented for each sequence as binary characters
           with values of '0', '1', or '-' (for missing data)
           Sequences that have no constraint may have a NULL string
        */
        int nConstraints;
        char **constraintSeqs;

        /* The profile data structures */
        int maxnode;			/* The next index to allocate */
        int maxnodes;			/* Space allocated in data structures below */
        //profile_t **profiles;         /* Profiles of leaves and intermediate nodes */
        numeric_t *diameter;		/* To correct for distance "up" from children (if any) */
        numeric_t *varDiameter;		/* To correct variances for distance "up" */
        numeric_t *selfdist;		/* Saved for use in some formulas */
        numeric_t *selfweight;		/* Saved for use in some formulas */

        /* Average profile of all active nodes, the "outprofile"
         * If all inputs are ungapped, this has weight 1 (not nSequences) at each position
         * The frequencies all sum to one (or that is implied by the eigen-representation)
         */
        //profile_t *outprofile;
        double totdiam;

        /* We sometimes use stale out-distances, so we remember what nActive was  */
        numeric_t *outDistances;		/* Sum of distances to other active (parent==-1) nodes */
        int *nOutDistActive;		/* What nActive was when this outDistance was computed */

        /* the inferred tree */
        int root;			/* index of the root. Unlike other internal nodes, it has 3 children */
        int *parent;			/* -1 or index of parent */
        //children_t *child;
        numeric_t *branchlength;		/* Distance to parent */
        numeric_t *support;		/* 1 for high-confidence nodes */

        /* auxilliary data for maximum likelihood (defaults to 1 category of rate=1.0) */
        //rates_t rates;
    };
}

#include "NeighbourJoining.tcc"

#endif
