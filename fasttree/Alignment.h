
#ifndef FASTTREE_ALIGNMENT_H
#define FASTTREE_ALIGNMENT_H

#include <vector>
#include <string>
#include "Options.h"

namespace fasttree {
    class Alignment {
    private:
        const Options &options;
        std::istream &fp;
        std::ostream &log;

    public:
        size_t nPos;
        std::vector<std::string> names;
        std::vector<std::string> seqs;

        Alignment(const Options &options, std::istream &fp, std::ostream &log);

        void readAlignment();

        void clearAlignmentSeqs();

        void clearAlignment();

    };

    /* Uniquify sequences in an alignment -- map from indices
       in the alignment to unique indicies in a NJ_t
    */
    struct Uniquify {
        std::vector<size_t> uniqueFirst;       /* iUnique -> iAln */
        std::vector<int64_t> alnNext;          /* iAln -> next, or -1  */
        std::vector<int64_t> alnToUniq;        /* iAln -> iUnique, or -1 if another was the exemplar */
        std::vector<std::string> uniqueSeq;    /* indexed by iUniq -- points to strings allocated elsewhere */

        Uniquify(const Alignment& aln);
    };
}

#endif
