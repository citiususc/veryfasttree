
#ifndef FASTTREE_NEIGHBOURJOINING_H
#define FASTTREE_NEIGHBOURJOINING_H

#include "DistanceMatrix.h"
#include "TransitionMatrix.h"
#include "HashTable.h"
#include "Alignment.h"
#include <vector>
#include <iostream>


namespace fasttree {

    template<typename Precision, template<typename> typename Operations>
    class NeighbourJoining {
    public:
        NeighbourJoining(Options &options, std::ostream &log, std::vector<std::string> &seqs, size_t nPos,
                         std::vector<std::string> &constraintSeqs,
                         DistanceMatrix<Precision> &distanceMatrix, TransitionMatrix<Precision> &transmat);

        void printDistances(std::vector<std::string> &names, std::ostream &out);

        /* ReadTree ignores non-unique leaves after the first instance.
           At the end, it prunes the tree to ignore empty children and it
           unroots the tree if necessary.
        */
        void readTree(Uniquify &unique, HashTable &hashnames, std::istream &fpInTree);

        void printNJ(std::ostream &out, std::vector<std::string> &names, Uniquify &unique, bool bShowSupport);

        /* Searches the visible set */
        void fastNJ();

    private:
        typedef Precision numeric_t;

        struct Profile {
            /* alignment profile */
            std::vector<numeric_t> weights;
            std::string codes;
            std::vector<numeric_t> vectors;        /* empty if no non-constant positions, e.g. for leaves */
            std::vector<numeric_t> codeDist;        /* Optional -- distance to each code at each position */

            /* constraint profile */
            std::vector<size_t> nOn;
            std::vector<size_t> nOff;

            numeric_t nGaps; /*precalculated in construction*/

            Profile(size_t nPos, size_t nConstraints);

            Profile();
        };

        struct Rates {
            std::vector<numeric_t> rates;    /* 1 per rate category */
            std::vector<size_t> ratecat;    /* 1 category per position */

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

        /* A visible node is a pair of nodes i, j such that j is the best hit of i,
           using the neighbor-joining criterion, at the time the comparison was made,
           or approximately so since then.

           Note that variance = dist because in BIONJ, constant factors of variance do not matter,
           and because we weight ungapped sequences higher naturally when averaging profiles,
           so we do not take this into account in the computation of "lambda" for BIONJ.

           For the top-hit list heuristic, if the top hit list becomes "too short",
           we store invalid entries with i=j=-1 and dist/criterion very high.
        */
        struct Besthit {
            int64_t i, j;
            numeric_t weight;        /* Total product of weights (maximum value is nPos)
                                           This is needed for weighted joins and for pseudocounts,
                                           but not in most other places.
                                           For example, it is not maintained by the top hits code */
            numeric_t dist;            /* The uncorrected distance (includes diameter correction) */
            numeric_t criterion;    /* changes when we update the out-profile or change nActive */
        };


        Operations<Precision> operations;
        std::ostream &log;
        Options &options;
        /* The input */
        size_t nPos;
        std::vector<std::string> &seqs;    /* the aligment sequences array (not reallocated) */
        DistanceMatrix<Precision> &distanceMatrix; /* a ref, not set if using %identity distance */
        TransitionMatrix<Precision> &transmat; /* a ref, not set for Jukes-Cantor */
        /* Topological constraints are represented for each sequence as binary characters
           with values of '0', '1', or '-' (for missing data)
           Sequences that have no constraint may have a NULL string
        */
        std::vector<std::string> &constraintSeqs;

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

        double logCorrect(double dist);

        void seqsToProfiles();

        /* outProfile() computes the out-profile */
        void outProfile(Profile &out);

        /* Use out-profile and NJ->totdiam to recompute out-distance for node iNode
           Only does this computation if the out-distance is "stale" (nOutDistActive[iNode] != nActive)
           Note "IN/UPDATE" for NJ always means that we may update out-distances but otherwise
           make no changes.
        */
        void setOutDistance(size_t iNode, size_t nActive);

        /* only handles leaf sequences */
        size_t nGaps(size_t i);

        /* E.g. GET_FREQ(profile,iPos,iVector)
           Gets the next element of the vectors (and updates iVector), or
           returns NULL if we didn't store a vector
        */
        inline numeric_t *getFreq(Profile &p, size_t &i, size_t &ivector);

        /* Adds (or subtracts, if weight is negative) fIn/codeIn from fOut
           fOut is assumed to exist (as from an outprofile)
           do not call unless weight of input profile > 0
         */
        void addToFreq(numeric_t fOut[], double weight, size_t codeIn, numeric_t fIn[]);

        /* Divide the vector (of length nCodes) by a constant
           so that the total (unrotated) frequency is 1.0
        */
        void normalizeFreq(numeric_t freq[]);

        /* Allocate, if necessary, and recompute the codeDist*/
        void setCodeDist(Profile &profile);

        /* f1 can be NULL if code1 != NOCODE, and similarly for f2
           Or, if (say) weight1 was 0, then can have code1==NOCODE *and* f1==NULL
           In that case, returns an arbitrary large number.
        */
        double profileDistPiece(size_t code1, size_t code2, numeric_t f1[], numeric_t f2[], numeric_t codeDist2[]);

        /* ProfileDist and SeqDist only set the dist and weight fields
           If using an outprofile, use the second argument of ProfileDist
           for better performance.

           These produce uncorrected distances.
        */
        void profileDist(Profile &profile1, Profile &profile2, Besthit &hit);

        void seqDist(std::string &codes1, std::string &codes2, Besthit &hit);


        /* Set a node's profile from its children.
           Deletes the previous profile if it exists
           Use -1.0 for a balanced join
           Fails unless the node has two children (e.g., no leaves or root)
        */
        void setProfile(size_t node, double weight1);

        /* AverageProfile is used to do a weighted combination of nodes
           when doing a join. If weight is negative, then the value is ignored and the profiles
           are averaged. The weight is *not* adjusted for the gap content of the nodes.
           Also, the weight does not affect the representation of the constraints
        */
        void averageProfile(Profile &out, Profile &profile1, Profile &profile2, double weight1);

        void readTreeRemove(std::vector<int64_t>& parents, std::vector<Children>& children, size_t node);

        /* A token is one of ():;, or an alphanumeric string without whitespace
            Any whitespace between tokens is ignored */
        bool readTreeToken(std::istream &fpInTree, std::string &buf);

        void
        readTreeAddChild(size_t parent, size_t child, std::vector<int64_t> &parents, std::vector<Children> &children);

        void readTreeMaybeAddLeaf(size_t parent, std::string &name, HashTable &hashnames, Uniquify &unique,
                                  std::vector<int64_t> &parents, std::vector<Children> &children);

        void readTreeError(const std::string &err, const std::string &token);

        /* returns new node, or -1 if nothing left to do. Use root for the first call.
           Will return every node and then root.
           Uses postorder tree traversal (depth-first search going down to leaves first)
           Keeps track of which nodes are visited, so even after an NNI that swaps a
           visited child with an unvisited uncle, the next call will visit the
           was-uncle-now-child. (However, after SPR moves, there is no such guarantee.)

           If pUp is not NULL, then, if going "back up" through a previously visited node
           (presumably due to an NNI), then it will return the node another time,
           with *pUp = true.
        */
        size_t TraversePostorder(int64_t lastnode, std::vector<bool> &traversal, bool *pUp);
    };
}

#include "NeighbourJoining.tcc"

#endif
