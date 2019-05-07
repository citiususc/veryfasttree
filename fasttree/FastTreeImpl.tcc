
#ifndef FASTTREE_FASTTREEIMPL_TCC
#define FASTTREE_FASTTREEIMPL_TCC

#include "FastTreeImpl.h"
#include <sstream>
#include "NeighbourJoining.h"
#include "assert.h"

#define AbsFastTreeImpl(...) \
template<typename Precision, template<class> class Operations> \
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
    for (int64_t iAln = 0; iAln < options.nAlign; iAln++) {
        aln.readAlignment();

        if (aln.seqs.empty()) {
            throw std::invalid_argument("No alignment sequences");
        }

        if (!options.logFileName.empty()) {
            log << strformat("Read %d sequences, %d positions", aln.seqs.size(), aln.nPos) << std::endl;
        }

        progressReport.print("Read alignment");


        /* Check that all names in alignment are unique */
        HashTable hashnames(aln.names, true);

        /* Make a list of unique sequences -- note some lists are bigger than required */
        progressReport.print("Hashed the names");
        if (options.makeMatrix) {
            std::vector<std::string> constraintSeqs;
            NeighbourJoining <Precision, Operations> nj(options, log, progressReport, aln.seqs, aln.nPos,
                                                        constraintSeqs, distanceMatrix, transmat);
            nj.printDistances(aln.names, output);
        } else {
            /* reset counters*/
            options.debug.reset();

            Uniquify unique(aln);
            progressReport.print("Identified unique sequences");

            /* read constraints */
            std::vector<std::string> uniqConstraints;
            if (fpConstraints) {
                Alignment constraints(options, fpConstraints, log);
                constraints.readAlignment();
                if (constraints.seqs.size() < 4) {
                    log << "Warning: constraints file with less than 4 sequences ignored:" << std::endl;
                    log << strformat("alignment #%d in ", iAln + 1) << options.constraintsFile << std::endl;
                } else {
                    alnToConstraints(uniqConstraints, constraints, unique, hashnames);
                    progressReport.print("Read the constraints");
                }
            }    /* end load constraints */

            if (options.nCodes == 20) {
                if (options.bUseLg) {
                    transmat.createTransitionMatrixLG08(options);
                } else if (options.bUseWag) {
                    transmat.createTransitionMatrixWAG01(options);
                } else {
                    transmat.createTransitionMatrixJTT92(options);
                }
            } else if (options.nCodes == 4 && options.bUseGtr && (options.bUseGtrRates || options.bUseGtrFreq)) {
                transmat.createGTR(options, options.gtrrates, options.gtrfreq);
            }
            NeighbourJoining <Precision, Operations> nj(options, log, progressReport, unique.uniqueSeq, aln.nPos,
                                                        uniqConstraints, distanceMatrix, transmat);
            if (options.verbose > 2) {
                log << strformat("read %s seqs %d (%d unique) positions %d nameLast %s seqLast %s",
                                 options.inFileName.empty() ? "standard input" : options.inFileName.c_str(),
                                 aln.seqs.size(), unique.uniqueSeq.size(),
                                 aln.nPos, aln.names[aln.seqs.size() - 1], aln.seqs[aln.seqs.size() - 1]) << std::endl;
            }
            aln.clearAlignmentSeqs(); /*no longer needed*/
            if (fpInTree) {
                if (options.intree1) {
                    fpInTree.clear();
                    fpInTree.seekg(0, std::ios::beg);
                }
                nj.readTree(unique, hashnames, fpInTree);
                if (options.verbose > 2) {
                    log << "Read tree from " << options.intreeFile << std::endl;
                }
                if (options.verbose > 2) {
                    nj.printNJ(log, aln.names, unique, false);
                }
            } else {
                nj.fastNJ();
            }
            nj.logTree("NJ", 0, aln.names, unique);

            int64_t nniToDo = options.nni == -1 ? (int64_t) (0.5 + 4.0 * std::log(aln.seqs.size()) / std::log(2)) :
                              options.nni;
            int64_t sprRemaining = options.spr;
            int64_t MLnniToDo = (options.MLnni != -1) ? options.MLnni :
                                (int64_t) (0.5 + 2.0 * std::log(aln.seqs.size()) / std::log(2));
            if (options.verbose > 0) {
                if (!fpInTree) {
                    log << strformat("Initial topology in %.2f seconds", progressReport.clockDiff()) << std::endl;
                }
                if (options.spr > 0 || nniToDo > 0 || MLnniToDo > 0) {
                    log << strformat("Refining topology: %d rounds ME-NNIs, %d rounds ME-SPRs, %d rounds ML-NNIs",
                                     nniToDo, options.spr, MLnniToDo) << std::endl;
                }
            }


            //TODO implement
        }
    }

}

AbsFastTreeImpl(void)::alnToConstraints(std::vector<std::string> &uniqConstraints, Alignment &constraints,
                                        Uniquify &unique, HashTable &hashnames) {
    /* look up constraints as names and map to unique-space */
    uniqConstraints.reserve(unique.uniqueSeq.size());

    for (int64_t i = 0; i < (int64_t) constraints.seqs.size(); i++) {
        auto &name = constraints.names[i];
        auto &constraintSeq = constraints.seqs[i];
        auto hi = hashnames.find(name);
        if (hi == nullptr) {
            throw std::invalid_argument(
                    strformat("Sequence %s from constraints file is not in the alignment", name.c_str()));
        }
        int64_t iSeqNonunique = *hi;
        assert(iSeqNonunique >= 0 && iSeqNonunique < (int64_t) unique.alnToUniq.size());
        int64_t iSeqUnique = unique.alnToUniq[iSeqNonunique];
        assert(iSeqUnique >= 0 && iSeqUnique < (int64_t) unique.uniqueSeq.size());
        if (!uniqConstraints[iSeqUnique].empty()) {
            /* Already set a constraint for this group of sequences!
            Warn that we are ignoring this one unless the constraints match */
            if (uniqConstraints[iSeqUnique] != constraintSeq) {
                log << strformat("Warning: ignoring constraints for %s:", name.c_str()) << std::endl;
                log << constraintSeq << std::endl;
                log << "Another sequence has the same sequence but different constraints" << std::endl;
            }
        } else {
            uniqConstraints[iSeqUnique] = constraintSeq;
        }
    }
}


#endif