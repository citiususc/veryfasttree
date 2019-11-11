
#ifndef VERYFASTTREE_VERYFASTTREE_H
#define VERYFASTTREE_VERYFASTTREE_H

#include "VeyFastTreeImpl.h"
#include "tsl/robin_map.h"

namespace veryfasttree {
    class VeryFastTree {
    private:
        Options options;

        void settings(std::ostream &log);
        void configOpenMP();

    public:
        VeryFastTree(const Options &options);

        void run(std::istream &input, std::ostream &output, std::ostream &log);

    };
}

#endif
