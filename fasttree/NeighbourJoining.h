
#ifndef FASTTREE_NEIGHBOURJOINING_H
#define FASTTREE_NEIGHBOURJOINING_H

#include "DistanceMatrix.h"
#include "TransitionMatrix.h"
#include "HashTable.h"
#include "Alignment.h"
#include <vector>
#include <iostream>


namespace fasttree {

    template<typename Precision, template<class> class Operations>
    class NeighbourJoining {
    public:
        NeighbourJoining(Options &options, std::ostream &log, ProgressReport &progressReport,
                         std::vector<std::string> &seqs, int64_t nPos,
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

        void logMLRates();

        void logTree(const std::string &format, int64_t i, std::vector<std::string> &names, Uniquify &unique);

        /* nni_stats_t is meaningless for leaves and root, so all of those entries
           will just be high (for age) or 0 (for delta)
        */
        struct NNIStats {
            int64_t age;        /* number of rounds since this node was modified by an NNI */
            int64_t subtreeAge; /* number of rounds since self or descendent had a significant improvement */
            double delta;        /* improvement in score for this node (or 0 if no change) */
            double support;     /* improvement of score for self over better of alternatives */
        };

        struct SplitCount {
            int64_t nBadSplits;
            int64_t nConstraintViolations;
            int64_t nBadBoth;
            int64_t nSplits;
            /* How much length would be reduce or likelihood would be increased by the
               best NNI we find (the worst "miss") */
            double dWorstDeltaUnconstrained;
            double dWorstDeltaConstrained;
        };

        /* One round of nearest-neighbor interchanges according to the
           minimum-evolution or approximate maximum-likelihood criterion.
           If doing maximum likelihood then this modifies the branch lengths.
           age is the # of rounds since a node was NNId
           Returns the # of topological changes performed
        */
        int64_t
        NNI(int64_t iRound, int64_t nRounds, bool useML, std::vector<NNIStats> &stats, double &maxDeltaCriterion);

        /* Recomputes all branch lengths by minimum evolution criterion*/
        void updateBranchLengths();

        /* Recomputes all branch lengths and, optionally, internal profiles */
        double treeLength(bool recomputeProfiles);

        void initNNIStats(std::vector<NNIStats> &stats);

        /* One round of subtree-prune-regraft moves (minimum evolution) */
        void SPR(int64_t maxSPRLength, int64_t iRound, int64_t nRounds);

    private:
        typedef Precision numeric_t;
        typedef Operations<Precision> op_t;

        struct Profile {
            /* alignment profile */
            std::vector<numeric_t> weights;
            std::string codes;
            std::vector<numeric_t> vectors;        /* empty if no non-constant positions, e.g. for leaves */
            std::vector<numeric_t> codeDist;        /* Optional -- distance to each code at each position */

            /* constraint profile */
            std::vector<int64_t> nOn;
            std::vector<int64_t> nOff;

            numeric_t nGaps; /*precalculated in construction*/

            Profile(int64_t nPos, int64_t nConstraints);

            Profile();
        };

        struct Rates {
            std::vector<numeric_t> rates;    /* 1 per rate category */
            std::vector<int64_t> ratecat;    /* 1 category per position */

            /* Allocate or reallocate the rate categories, and set every position
               to category 0 and every category's rate to 1.0
               If nRateCategories=0, just deallocate
            */
            Rates(int64_t nRateCategories, int64_t nPos);
        };

        struct Children {
            int nChild = 0;
            int64_t child[3];
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


        /* Keep track of hits for the top-hits heuristic without wasting memory
           j = -1 means empty
           If j is an inactive node, this may be replaced by that node's parent (and dist recomputed)
         */
        struct Hit {
            int64_t j;
            numeric_t dist;
        };

        struct TopHitsList {
            std::vector<Hit> hits;  /* the allocated and desired size; some of them may be empty */
            int64_t hitSource;        /* where to refresh hits from if a 2nd-level top-hit list, or -1 */
            int64_t age;                /* number of joins since a refresh */
        };

        struct TopHits {
            int64_t m;             /* size of a full top hits list, usually sqrt(N) */
            int64_t q;             /* size of a 2nd-level top hits, usually sqrt(m) */
            int64_t maxnodes;
            std::vector<TopHitsList> topHitsLists; /* one per node */
            std::vector<Hit> visible;        /* the "visible" (very best) hit for each node */

            /* The top-visible set is a subset, usually of size m, of the visible set --
               it is the set of joins to select from
               Each entry is either a node whose visible set entry has a good (low) criterion,
               or -1 for empty, or is an obsolete node (which is effectively the same).
               Whenever we update the visible set, should also call UpdateTopVisible()
               which ensures that none of the topvisible set are stale (that is, they
               all point to an active node).
            */
            std::vector<int64_t> topvisible; /* nTopVisible = m * topvisibleMult */

            int64_t topvisibleAge;        /* joins since the top-visible list was recomputed */


            /* 1 lock to read or write any top hits list, no thread grabs more than one */
            //omp_lock_t *locks;

            TopHits(const Options &options, int64_t maxnodes, int64_t m);

            TopHits();

        };

        /* Describes which switch to do */
        enum NNI {
            ABvsCD, ACvsBD, ADvsBC
        };


        Operations<Precision> operations;
        std::ostream &log;
        Options &options;
        ProgressReport &progressReport;
        /* The input */
        int64_t nPos;
        std::vector<std::string> &seqs;    /* the aligment sequences array (not reallocated) */
        DistanceMatrix<Precision> &distanceMatrix; /* a ref, not set if using %identity distance */
        TransitionMatrix<Precision> &transmat; /* a ref, not set for Jukes-Cantor */
        /* Topological constraints are represented for each sequence as binary characters
           with values of '0', '1', or '-' (for missing data)
           Sequences that have no constraint may have a NULL string
        */
        std::vector<std::string> &constraintSeqs;

        /* The profile data structures */
        int64_t maxnode;            /* The next index to allocate */
        int64_t maxnodes;            /* Space allocated in data structures below */
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
        std::vector<int64_t> nOutDistActive;        /* What nActive was when this outDistance was computed */

        /* the inferred tree */
        int64_t root;                   /* index of the root. Unlike other internal nodes, it has 3 children */
        std::vector<int64_t> parent;    /* -1 or index of parent */
        std::vector<Children> child;
        std::vector<numeric_t> branchlength;        /* Distance to parent */
        std::vector<numeric_t> support;        /* 1 for high-confidence nodes */

        /* auxilliary data for maximum likelihood (defaults to 1 category of rate=1.0) */
        Rates rates;

        double logCorrect(double dist);

        /* Print topology using node indices as node names */
        void printNJInternal(std::ostream &out, bool useLen);

        void seqsToProfiles();

        /* This recomputes the criterion, or returns false if the visible node
           is no longer active.
        */
        bool getVisible(int64_t nActive, TopHits &tophits, int64_t iNode, Besthit &visible);

        int64_t activeAncestor(int64_t iNode);

        /* Compute the constraint penalty for a join. This is added to the "distance"
           by SetCriterion */
        int64_t joinConstraintPenalty(int64_t node1, int64_t node2);

        int64_t joinConstraintPenaltyPiece(int64_t node1, int64_t node2, int64_t iConstraint);

        /* Helper function for computing the number of constraints violated by
           a split, represented as counts of on and off on each side */
        int64_t splitConstraintPenalty(int64_t nOn1, int64_t nOff1, int64_t nOn2, int64_t nOff2);

        /* outProfile() computes the out-profile,
         * Profile_t can be Profile or Profile*
         * */
        template<typename Profile_t>
        void outProfile(Profile &out, std::vector<Profile_t> &_profiles);

        void updateOutProfile(Profile &out, Profile &old1, Profile &old2, Profile &_new, int nActiveOld);

        /* Use out-profile and NJ->totdiam to recompute out-distance for node iNode
           Only does this computation if the out-distance is "stale" (nOutDistActive[iNode] != nActive)
           Note "IN/UPDATE" for NJ always means that we may update out-distances but otherwise
           make no changes.
        */
        void setOutDistance(int64_t iNode, int64_t nActive);

        /* Always sets join->criterion; may update NJ->outDistance and NJ->nOutDistActive,
           assumes join's weight and distance are already set,
           and that the constraint penalty (if any) is included in the distance
        */
        void setCriterion(int64_t nActive, Besthit &join);

        /* Computes weight and distance (which includes the constraint penalty)
           and then sets the criterion (maybe update out-distances)
        */
        void setDistCriterion(int64_t nActive, Besthit &join);

        /* If join->i or join->j are inactive nodes, replaces them with their active ancestors.
           After doing this, if i == j, or either is -1, sets weight to 0 and dist and criterion to 1e20
              and returns false (not a valid join)
           Otherwise, if i or j changed, recomputes the distance and criterion.
           Note that if i and j are unchanged then the criterion could be stale
           If bUpdateDist is false, and i or j change, then it just sets dist to a negative number
        */
        bool updateBestHit(int64_t nActive, Besthit &join, bool bUpdateDist);

        /* Returns the resulting log likelihood. Optionally returns whether other
           topologies should be abandoned, based on the difference between AB|CD and
           the "star topology" (AB|CD with a branch length of MLMinBranchLength) exceeding
           closeLogLkLimit.
           If bStarTest is passed in, it only optimized the internal branch if
           the star test is true. Otherwise, it optimized all 5 branch lengths
           in turn.
         */
        double MLQuartetOptimize(Profile &pA, Profile &pB, Profile &pC, Profile &pD, double branch_lengths[5],
                                 bool *pStarTest, double *site_likelihoods);

        /* Update the profile of node and its ancestor, and delete nearby out-profiles */
        void updateForNNI(int64_t node, std::unique_ptr<Profile> upProfiles[], bool useML);

        /* Sets NJ->parent[newchild] and replaces oldchild with newchild
           in the list of children of parent
        */
        void replaceChild(int64_t parent, int64_t oldchild, int64_t newchild);

        /* only handles leaf sequences */
        int64_t nGaps(int64_t i);

        /* node is the parent of AB, sibling of C
           node cannot be root or a leaf
           If node is the child of root, then D is the other sibling of node,
           and the 4th profile is D's profile.
           Otherwise, D is the parent of node, and we use its upprofile
           Call this with profiles=NULL to get the nodes, without fetching or
           computing profiles
        */
        void setupABCD(int64_t node, Profile *profiles[4], std::unique_ptr<Profile> upProfiles[], int64_t nodeABCD[4],
                       bool useML);

        int64_t sibling(int64_t node); /* At root, no unique sibling so returns -1 */
        void rootSiblings(int64_t node, /*OUT*/int64_t sibs[2]);

        /* JC probability of nucleotide not changing, for each rate category */
        void pSameVector(double length, std::vector<double> &pSame);

        /* JC probability of nucleotide not changing, for each rate category */
        void pDiffVector(std::vector<double> &pSame, std::vector<double> &pDiff);

        /* expeigen[iRate*nCodes + j] = exp(length * rate iRate * eigenvalue j) */
        void expEigenRates(double length, std::vector<numeric_t> &expeigenRates);

        /* E.g. GET_FREQ(profile,iPos,iVector)
           Gets the next element of the vectors (and updates iVector), or
           returns NULL if we didn't store a vector
        */
        inline numeric_t *getFreq(Profile &p, int64_t &i, int64_t &ivector);

        /* Adds (or subtracts, if weight is negative) fIn/codeIn from fOut
           fOut is assumed to exist (as from an outprofile)
           do not call unless weight of input profile > 0
         */
        void addToFreq(numeric_t fOut[], double weight, int64_t codeIn, numeric_t fIn[]);

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
        double profileDistPiece(int64_t code1, int64_t code2, numeric_t f1[], numeric_t f2[], numeric_t codeDist2[]);

        /* ProfileDist and SeqDist only set the dist and weight fields
           If using an outprofile, use the second argument of ProfileDist
           for better performance.

           These produce uncorrected distances.
        */
        void profileDist(Profile &profile1, Profile &profile2, Besthit &hit);

        /* P(A & B | len) = P(B | A, len) * P(A)
           If site_likelihoods is present, multiplies those values by the site likelihood at each point
           (Note it does not handle underflow)
         */
        double pairLogLk(Profile &p1, Profile &p2, double length, double site_likelihoods[]);

        /* Computes all pairs of profile distances, applies pseudocounts
           if pseudoWeight > 0, and applies log-correction if logdist is true.
           The lower index is compared to the higher index, e.g. for profiles
           A,B,C,D the comparison will be as in quartet_pair_t
        */
        enum QuartetPair {
            qAB, qAC, qAD, qBC, qBD, qCD
        };

        struct QuartetOpt {
            int64_t nEval;            /* number of likelihood evaluations */
            /* The pair to optimize */
            Profile *pair1;
            Profile *pair2;
        };

        double pairNegLogLk(double x, QuartetOpt &qo);

        /* Branch lengths for 4-taxon tree ((A,B),C,D); I means internal */
        enum QuartetLength {
            LEN_A, LEN_B, LEN_C, LEN_D, LEN_I
        };

        void correctedPairDistances(Profile *profiles[], int nProfiles, double distances[6]);

        /* output is indexed by nni_t
           To ensure good behavior while evaluating a subtree-prune-regraft move as a series
           of nearest-neighbor interchanges, this uses a distance-ish model of constraints,
           as given by PairConstraintDistance(), rather than
           counting the number of violated splits (which is what FastTree does
           during neighbor-joining).
           Thus, penalty values may well be >0 even if no constraints are violated, but the
           relative scores for the three NNIs will be correct.
         */
        void quartetConstraintPenalties(Profile *profiles[4], double penalty[3]);

        double pairConstraintDistance(int64_t nOn1, int64_t nOff1, int64_t nOn2, int64_t nOff2);

        /* the split is consistent with the constraint if any of the profiles have no data
           or if three of the profiles have the same uniform value (all on or all off)
           or if AB|CD = 00|11 or 11|00 (all uniform)
         */
        bool splitViolatesConstraint(Profile *profiles[4], int64_t iConstraint);

        /* If false, no values were set because this constraint was not relevant.
           output is for the 3 splits
        */
        bool quartetConstraintPenaltiesPiece(Profile *profiles[4], int64_t iConstraint, double piece[3]);

        void seqDist(std::string &codes1, std::string &codes2, Besthit &hit);


        /* Set a node's profile from its children.
           Deletes the previous profile if it exists
           Use -1.0 for a balanced join
           Fails unless the node has two children (e.g., no leaves or root)
        */
        void setProfile(int64_t node, double weight1);

        /* AverageProfile is used to do a weighted combination of nodes
           when doing a join. If weight is negative, then the value is ignored and the profiles
           are averaged. The weight is *not* adjusted for the gap content of the nodes.
           Also, the weight does not affect the representation of the constraints
        */
        void averageProfile(Profile &out, Profile &profile1, Profile &profile2, double weight1);

        /* PosteriorProfile() is like AverageProfile() but it computes posterior probabilities
           rather than an average
        */
        void posteriorProfile(Profile &out, Profile &profile1, Profile &profile2, double len1, double len2);

        void readTreeRemove(std::vector<int64_t> &parents, std::vector<Children> &children, int64_t node);

        /* A token is one of ():;, or an alphanumeric string without whitespace
            Any whitespace between tokens is ignored */
        bool readTreeToken(std::istream &fpInTree, std::string &buf);

        void
        readTreeAddChild(int64_t parent, int64_t child, std::vector<int64_t> &parents, std::vector<Children> &children);

        void readTreeMaybeAddLeaf(int64_t parent, std::string &name, HashTable &hashnames, Uniquify &unique,
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
        int64_t traversePostorder(int64_t lastnode, std::vector<bool> &traversal, bool *pUp);

        Profile *getUpProfile(std::unique_ptr<Profile> upProfiles[], int64_t outnode, bool useML);

        /* Recomputes the profile for a node, presumably to reflect topology changes
           If bionj is set, does a weighted join -- which requires using upProfiles
           If useML is set, computes the posterior probability instead of averaging
         */
        void recomputeProfile(std::unique_ptr<Profile> upProfiles[], int64_t node, int64_t useML);

        /* Recompute profiles going up from the leaves, using the provided distance matrix
           and unweighted joins
        */
        void recomputeProfiles();

        void recomputeMLProfiles();

        /* If bionj is set, computes the weight to be given to A when computing the
           profile for the ancestor of A and B. C and D are the other profiles in the quartet
           If bionj is not set, returns -1 (which means unweighted in AverageProfile).
           (A and B are the first two profiles in the array)
        */
        double quartetWeight(Profile *profiles[4]);

        /* Returns a list of nodes, starting with node and ending with root */
        void pathToRoot(int64_t node, std::vector<int64_t> &path);

        /* The allhits list contains the distances of the node to all other active nodes
           This is useful for the "reset" improvement to the visible set
           Note that the following routines do not handle the tophits heuristic
           and assume that out-distances are up to date.
        */
        void setBestHit(int64_t node, int64_t nActive, Besthit &bestjoin, Besthit allhits[]);

        void exhaustiveNJSearch(int64_t nActive, Besthit &bestjoin);

        /* Searches the visible set */
        void fastNJSearch(int64_t nActive, std::vector<Besthit> &visible, Besthit &bestjoin);

        /* Subroutines for handling the tophits heuristic */

        /* Before we do any joins -- sets tophits and visible
           NJ may be modified by setting out-distances
         */
        void setAllLeafTopHits(TopHits &tophits);

        /* Find the best join to do. */
        void topHitNJSearch(int64_t nActive, TopHits &tophits, Besthit &bestjoin);

        /* Returns the best hit within top hits
           NJ may be modified because it updates out-distances if they are too stale
           Does *not* update visible set
        */
        void getBestFromTopHits(int64_t iNode, int64_t nActive, TopHits &tophits, Besthit &bestjoin);

        /* visible set is modifiable so that we can reset it more globally when we do
           a "refresh", but we also set the visible set for newnode and do any
           "reset" updates too. And, we update many outdistances.
         */
        void topHitJoin(int64_t newnode, int64_t nActive, TopHits &tophits);

        /* Sort the input besthits by criterion
           and save the best nOut hits as a new array in top_hits_lists
           Does not update criterion or out-distances
           Ignores (silently removes) hit to self
           Saved list may be shorter than requested if there are insufficient entries
        */
        void sortSaveBestHits(int64_t iNode, std::vector<Besthit> &besthits, int64_t nIn, int64_t nOut,
                              TopHits &tophits);

        /* Given candidate hits from one node, "transfer" them to another node:
           Stores them in a new place in the same order
           searches up to active nodes if hits involve non-active nodes
           If update flag is set, it also recomputes distance and criterion
           (and ensures that out-distances are updated); otherwise
           it sets dist to -1e20 and criterion to 1e20

         */
        void transferBestHits(int64_t nActive, int64_t iNode, std::vector<Besthit> &oldhits, int64_t nOldHits,
                              Besthit newhits[], bool updateDistance);

        /* Create best hit objects from 1 or more hits. Do not update out-distances or set criteria */
        void hitsToBestHits(std::vector<Hit> &hits, int64_t iNode, Besthit newhits[]);

        void hitToBestHit(int64_t i, Hit &hit, Besthit &out);

        /* Given a set of besthit entries,
           look for improvements to the visible set of the j entries.
           Updates out-distances as it goes.
           Also replaces stale nodes with this node, because a join is usually
           how this happens (i.e. it does not need to walk up to ancestors).
           Note this calls UpdateTopVisible() on any change
        */
        void updateVisible(int64_t nActive, std::vector<Besthit> &tophitsNode, TopHits &tophits);

        /* Update the top-visible list to perhaps include this hit (O(sqrt(N)) time) */
        void updateTopVisible(int64_t nActive, int64_t iNode, Hit &hit, TopHits &tophits);

        /* Recompute the top-visible subset of the visible set */
        void resetTopVisible(int64_t nActive, TopHits &tophits);

        /* Make a shorter list with only unique entries.
           Replaces any "dead" hits to nodes that have parents with their active ancestors
           and ignores any that become dead.
           Updates all criteria.
           Combined gets sorted by i & j
           The returned list is allocated to nCombined even though only *nUniqueOut entries are filled
        */
        void uniqueBestHits(int64_t nActive, std::vector<Besthit> &combined, std::vector<Besthit> &out);

        /* The three internal branch lengths or log likelihoods*/
        enum NNI chooseNNI(Profile *profiles[4], double criteria[3]);

        /* length[] is ordered as described by quartet_length_t, but after we do the swap
           of B with C (to give AC|BD) or B with D (to get AD|BC), if that is the returned choice
           bFast means do not consider NNIs if AB|CD is noticeably better than the star topology
           (as implemented by MLQuartetOptimize).
           If there are constraints, then the constraint penalty is included in criteria[]
        */
        enum NNI
        MLQuartetNNI(Profile *profiles[4], double criteria[3], /* The three potential quartet log-likelihoods */
                     numeric_t length[5], bool bFast);

        void optimizeAllBranchLengths();

        double treeLogLk(double site_loglk[]);

        double MLQuartetLogLk(Profile &pA, Profile &pB, Profile &pC, Profile &pD, double branchLengths[5],
                              double siteLikelihoods[]);

        /* Given a topology and branch lengths, estimate rates & recompute profiles */
        void setMLRates();

        /* Returns a set of nRateCategories potential rates; the caller must free it */
        void MLSiteRates(int64_t nRateCategories, std::vector<numeric_t> &rates);

        /* returns site_loglk so that
           site_loglk[nPos*iRate + j] is the log likelihood of site j with rate iRate
           The caller must free it.
        */
        void MLSiteLikelihoodsByRate(std::vector<numeric_t> &rates, std::vector<double> &site_loglk);

        /* One-dimensional minimization using brent's function, with
        a fractional and an absolute tolerance */
        template<typename Function, typename Data>
        double onedimenmin(double xmin, double xguess, double xmax, Function f, Data &data, double ftol, double atol,
                           double &fx, double &f2x);

        template<typename Function, typename Data>
        double brent(double ax, double bx, double cx, Function f, Data &data, double ftol, double atol,
                     double &foptx, double &f2optx, double fax, double fbx, double fcx);

        /************************Sort comparators************************/

        /* Helper function for sorting in SetAllLeafTopHits,
           and the global variables it needs
        */
        struct CompareSeeds {
            const std::vector<numeric_t> &outDistances;
            const std::vector<int64_t> &compareSeedGaps;

            CompareSeeds(const std::vector<numeric_t> &outDistances, const std::vector<int64_t> &compareSeedGaps);

            bool operator()(int64_t seed1, int64_t seed2) const;
        };

        struct CompareHitsByCriterion {
            bool operator()(const Besthit &hit1, const Besthit &hit2) const;
        };

        struct CompareHitsByIJ {
            bool operator()(const Besthit &hit1, const Besthit &hit2) const;
        };

    };
}

#include "NeighbourJoining.tcc"

#endif
