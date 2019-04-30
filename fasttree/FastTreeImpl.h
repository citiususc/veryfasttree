
#ifndef FASTTREE_FASTTREEIMPL_H
#define FASTTREE_FASTTREEIMPL_H


#include <memory>
#include <iostream>
#include <fstream>
#include "Options.h"
#include "Utils.h"
#include "DistanceMatrix.h"
#include "NeighbourJoining.h"

namespace fasttree {
    template<typename Precision, template<typename> typename Operations>
    class FastTreeImpl {
    private:
        Options &options;
        std::istream &input;
        std::ostream &output;
        std::ostream &log;

        std::ifstream fpConstraints;
        std::ifstream fpInTree;

        DistanceMatrix<Precision> distanceMatrix;
        TransitionMatrix<Precision> transmat;
        ProgressReport progressReport;

    public:
        FastTreeImpl(Options &options, std::istream &input, std::ostream &output, std::ostream &log);

        void run();
    };
}

#include "FastTreeImpl.tcc"

#endif
