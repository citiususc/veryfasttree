
#ifndef FASTTREE_FASTTREE_H
#define FASTTREE_FASTTREE_H

#include "FastTreeImpl.h"
#include "tsl/robin_map.h"

namespace fasttree {
    class FastTree {
    private:
        Options options;

        void prepare(std::istream &input, std::ostream &log);

    public:
        FastTree(const Options &options);

        void run(std::istream &input, std::ostream &output, std::ostream &log);

    };
}

#endif
