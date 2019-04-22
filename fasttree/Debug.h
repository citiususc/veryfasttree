
#ifndef FASTTREE_DEBUG_H
#define FASTTREE_DEBUG_H

#include <memory>

namespace fasttree {

    struct Debug {

        /* Performance and memory usage */
        size_t profileOps = 0;        /* Full profile-based distance operations */
        size_t outprofileOps = 0;        /* How many of profileOps are comparisons to outprofile */
        size_t seqOps = 0;        /* Faster leaf-based distance operations */
        size_t profileAvgOps = 0;        /* Number of profile-average steps */
        size_t nHillBetter = 0;        /* Number of hill-climbing steps */
        size_t nCloseUsed = 0;        /* Number of "close" neighbors we avoid full search for */
        size_t nClose2Used = 0;        /* Number of "close" neighbors we use 2nd-level top hits for */
        size_t nRefreshTopHits = 0;    /* Number of full-blown searches (interior nodes) */
        size_t nVisibleUpdate = 0;        /* Number of updates of the visible set */
        size_t nNNI = 0;            /* Number of NNI changes performed */
        size_t nSPR = 0;            /* Number of SPR changes performed */
        size_t nML_NNI = 0;        /* Number of max-lik. NNI changes performed */
        size_t nSuboptimalSplits = 0;    /* # of splits that are rejected given final tree (during bootstrap) */
        size_t nSuboptimalConstrained = 0; /* Bad splits that are due to constraints */
        size_t nConstraintViolations = 0;    /* Number of constraint violations */
        size_t nProfileFreqAlloc = 0;
        size_t nProfileFreqAvoid = 0;
        size_t szAllAlloc = 0;
        size_t mymallocUsed = 0;        /* useful allocations by mymalloc */
        size_t maxmallocHeap = 0;        /* Maximum of mi.arena+mi.hblkhd from mallinfo (actual mem usage) */
        size_t nLkCompute = 0;        /* # of likelihood computations for pairs of probability vectors */
        size_t nPosteriorCompute = 0;    /* # of computations of posterior probabilities */
        size_t nAAPosteriorExact = 0;    /* # of times compute exact AA posterior */
        size_t nAAPosteriorRough = 0;    /* # of times use rough approximation */
        size_t nStarTests = 0;        /* # of times we use star test to avoid testing an NNI */

        void reset() {
            (*this) = Debug();
        }
    };

    template<class Tp>
    class DebugAllocator {
    private :
        size_t maxHeapSize;
        size_t &maxHeapSizeRef;
        size_t heapSize;
        size_t &heapSizeRef;
        std::allocator<Tp> imp;

        template<typename>
        friend
        class DebugAllocator;

    public:

        typedef Tp value_type;
        typedef Tp *pointer;
        typedef const Tp *const_pointer;
        typedef Tp &reference;
        typedef const Tp &const_reference;
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;


        DebugAllocator() : maxHeapSizeRef(maxHeapSize), heapSizeRef(heapSize) {
            maxHeapSize = 0;
            heapSize = 0;
        }

        DebugAllocator(const DebugAllocator &alloc) :
                maxHeapSizeRef((size_t &) alloc.maxHeapSize),
                heapSizeRef((size_t &) alloc.heapSize) {}

        template<class U>
        DebugAllocator(const DebugAllocator<U> &alloc) :
                maxHeapSizeRef((size_t &) alloc.maxHeapSize),
                heapSizeRef((size_t &) alloc.heapSize) {}

        ~DebugAllocator() {
        }

        inline pointer allocate(size_t n, const void *hint = 0) {
            Tp *result = imp.allocate(n, hint);
            heapSizeRef += n;
            if (heapSizeRef > maxHeapSizeRef) {
                maxHeapSizeRef = heapSizeRef;
            }
            return result;
        }

        inline void deallocate(pointer p, std::size_t n) {
            imp.deallocate(p, n);
            heapSizeRef -= n;
        }

        inline size_t getMaxHeapSize() {
            return maxHeapSizeRef;
        }

        inline size_t getHeapSize() {
            return heapSizeRef;
        }

    };

}

#endif
