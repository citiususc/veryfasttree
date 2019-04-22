
#ifndef FASTTREE_FASTTREEIMPL_H
#define FASTTREE_FASTTREEIMPL_H

#include "Options.h"
#include <iostream>
#include <fstream>

namespace fasttree {
    template<typename Precision, template<typename> typename Operations, template<typename> typename Allocator>
    class FastTreeImpl {
    private:

        Options &options;
        std::ifstream fpIn;
        std::ifstream fpConstraints;
        std::ifstream fpInTree;


    public:
        FastTreeImpl(Options &options);

        void run(std::istream &input, std::ostream &output, std::ostream &log){}
    };
}

#include "FastTreeImpl.tcc"

#endif
