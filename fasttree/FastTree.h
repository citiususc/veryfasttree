
#ifndef FASTTREE_FASTTREE_H
#define FASTTREE_FASTTREE_H

#include "Options.h"
#include <iostream>
#include <fstream>

namespace fasttree {
    class FastTree {
    private:
        Options options;
        std::ifstream fpIn;
        std::ifstream fpConstraints;
        std::ifstream fpInTree;

        void prepare(std::istream &in, std::ostream &log);

    public:
        FastTree(const Options &options);

        void run(std::istream& input,  std::ostream& output, std::ostream &log);

    };
}

#endif
