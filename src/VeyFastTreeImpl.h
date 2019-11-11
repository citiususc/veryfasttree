
#ifndef VERYFASTTREE_VERYFASTTREEIMPL_H
#define VERYFASTTREE_VERYFASTTREEIMPL_H


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

namespace veryfasttree {
    template<typename Precision, template<class> class Operations>
    class VeyFastTreeImpl {
    public:
        VeyFastTreeImpl(Options &options, std::istream &input, std::ostream &output, std::ostream &log);

        void run();

    private:
        typedef Precision numeric_t;
        typedef Operations<Precision> op_t;

        Options &options;
        std::istream &input;
        std::ostream &output;
        std::ostream &log;

        std::ifstream fpConstraints;
        std::ifstream fpInTree;

        DistanceMatrix<Precision, op_t::ALIGNMENT> distanceMatrix;
        TransitionMatrix<Precision, op_t::ALIGNMENT> transmat;
        ProgressReport progressReport;

        /* Convert a constraint alignment to a list of sequences. The returned array is indexed
           by iUnique and points to values in the input alignment
        */
        void alnToConstraints(std::vector<std::string> &uniqConstraints, Alignment &constraints, Uniquify &unique,
                              HashTable &hashnames);

        void transMatToDistanceMat(DistanceMatrix<Precision, op_t::ALIGNMENT>& dmat);
    };
}

#include "VeryFastTreeImpl.tcc"

#endif
