
#ifndef FASTTREE_FASTTREE_H
#define FASTTREE_FASTTREE_H

#include "Options.h"
#include <iostream>

namespace fasttree {
    class FastTree {
    private:
        Options options;

    public:
        FastTree(const Options &options);

        void run(std::istream& input,  std::ostream& output, std::ostream &log);

    };
}

#endif
