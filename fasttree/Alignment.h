
#ifndef FASTTREE_ALIGNMENT_H
#define FASTTREE_ALIGNMENT_H

#include <vector>
#include <string>
#include "Options.h"

namespace fasttree {
    class Alignment {
    private:
        const Options& options;
        std::istream &fp;
        std::ostream &log;

        size_t _nPos;
        size_t _nSeq;
        std::vector<std::string> names;
        std::vector<std::string> seqs;
    public:
        Alignment(const Options &options, std::istream &fp, std::ostream &log);

        void readAlignment();

        void clearAlignmentSeqs();

        void clearAlignment();

        const size_t& nPos() const;

        const size_t& nSeq() const;

    };
}

#endif
