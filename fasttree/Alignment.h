
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
}

#endif
