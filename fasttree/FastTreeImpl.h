
#ifndef FASTTREE_FASTTREEIMPL_H
#define FASTTREE_FASTTREEIMPL_H


#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include "Options.h"
#include "Utils.h"
#include "Alignment.h"
#include "HashTable.h"
#include "DistanceMatrix.h"
#include "NeighbourJoining.h"

namespace fasttree {
    template<typename Precision, template<class> class Operations>
    class FastTreeImpl {
    public:
        FastTreeImpl(Options &options, std::istream &input, std::ostream &output, std::ostream &log);

        void run();

    private:
        typedef Precision numeric_t;

        Options &options;
        std::istream &input;
        std::ostream &output;
        std::ostream &log;

        std::ifstream fpConstraints;
        std::ifstream fpInTree;

        DistanceMatrix<Precision> distanceMatrix;
        TransitionMatrix<Precision> transmat;
        ProgressReport progressReport;

        /* Convert a constraint alignment to a list of sequences. The returned array is indexed
           by iUnique and points to values in the input alignment
        */
        void alnToConstraints(std::vector<std::string> &uniqConstraints, Alignment &constraints, Uniquify &unique,
                              HashTable &hashnames);

        void transMatToDistanceMat(DistanceMatrix<Precision>& dmat);
    };
}

#include "FastTreeImpl.tcc"

#endif
