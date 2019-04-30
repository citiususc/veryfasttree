
#ifndef FASTTREE_FASTTREEIMPL_TCC
#define FASTTREE_FASTTREEIMPL_TCC

#include <sstream>
#include <vector>
#include "FastTreeImpl.h"
#include "Alignment.h"
#include "NeighbourJoining.h"
#include "HashTable.h"
#include "assert.h"

#define AbsFastTreeImpl(...) \
template<typename Precision, template<typename> typename Operations> \
__VA_ARGS__ fasttree::FastTreeImpl<Precision, Operations>


AbsFastTreeImpl()::FastTreeImpl(Options &options, std::istream &input, std::ostream &output, std::ostream &log) :
        options(options), input(input), output(output), log(log), progressReport(options) {
    if (!options.matrixPrefix.empty()) {
        if (!options.useMatrix) {
            throw std::invalid_argument("Cannot use both -matrix and -nomatrix arguments!");
        }
        distanceMatrix.readDistanceMatrix(options, log);
        distanceMatrix.setupDistanceMatrix(options, log);
    } else if (options.useMatrix) {   // use default matrix
        assert(options.nCodes == 20);
        distanceMatrix.matrixBLOSUM45();
        distanceMatrix.setupDistanceMatrix(options, log);
    }

    if (!options.constraintsFile.empty()) {
        fpConstraints.open(options.constraintsFile);
        if (fpConstraints.fail()) {
            throw std::invalid_argument("Cannot read " + options.constraintsFile);
        }
    }

    if (!options.intreeFile.empty()) {
        fpInTree.open(options.intreeFile);
        if (fpInTree.fail()) {
            throw std::invalid_argument("Cannot read " + options.intreeFile);
        }
    }
}

AbsFastTreeImpl(void)::run() {
    Alignment aln(options, input, log);
    for (size_t iAln = 0; iAln < options.nAlign; iAln++) {
        aln.readAlignment();

        if (aln.seqs.empty()) {
            throw std::invalid_argument("No alignment sequences");
        }

        if (!options.logFileName.empty()) {
            log << strformat("Read %d sequences, %d positions", aln.seqs.size(), aln.nPos) << std::endl;
        }

        progressReport.print("Read alignment");


        /* Check that all names in alignment are unique */
        HashTable hashnames(aln.names);

        /* Make a list of unique sequences -- note some lists are bigger than required */
        progressReport.print("Hashed the names");
        if (options.makeMatrix) {
            std::vector<std::string> constraintSeqs;
            NeighbourJoining<Precision, Operations> nj(options, log, aln.seqs, aln.nPos, constraintSeqs, distanceMatrix,
                                                       transmat);
            nj.printDistances(aln.names, output);
        } else {
            //TODO implement
        }
    }

}


#endif