
#ifndef FASTTREE_FASTTREEIMPL_TCC
#define FASTTREE_FASTTREEIMPL_TCC

#include "FastTreeImpl.h"
#include <sstream>
#include "NeighbourJoining.h"
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
        HashTable hashnames(aln.names, true);

        /* Make a list of unique sequences -- note some lists are bigger than required */
        progressReport.print("Hashed the names");
        if (options.makeMatrix) {
            std::vector<std::string> constraintSeqs;
            NeighbourJoining <Precision, Operations> nj(options, log, aln.seqs, aln.nPos, constraintSeqs,
                                                        distanceMatrix,
                                                        transmat);
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
            NeighbourJoining <Precision, Operations> nj(options, log, unique.uniqueSeq, aln.nPos, uniqConstraints,
                                                        distanceMatrix,
                                                        transmat);
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
                //ReadTree(/*IN/OUT*/NJ, /*IN*/unique, /*IN*/hashnames, /*READ*/fpInTree);//TODO
                if (options.verbose > 2) {
                    log << "Read tree from " << options.intreeFile << std::endl;
                }
                if (options.verbose > 2) {
                    //PrintNJ(stderr, NJ, aln->names, unique, /*support*/false, bQuote);//TODO
                }
            } else {
                //FastNJ(NJ);//TODO
            }
            //LogTree("NJ", 0, fpLog, NJ, aln->names, unique, bQuote); //TODO

            /* profile-frequencies for the "up-profiles" in ReliabilityNJ take only diameter(Tree)*L*a
               space not N*L*a space, because we can free them as we go.
               And up-profile by their nature tend to be complicated.
               So save the profile-frequency memory allocation counters now to exclude later results.
            */

            //TODO implement
        }
    }

}

AbsFastTreeImpl(void)::alnToConstraints(std::vector<std::string> &uniqConstraints, Alignment &constraints,
                                        Uniquify &unique, HashTable &hashnames) {
    /* look up constraints as names and map to unique-space */
    uniqConstraints.reserve(unique.uniqueSeq.size());

    for (size_t i = 0; i < constraints.seqs.size(); i++) {
        auto &name = constraints.names[i];
        auto &constraintSeq = constraints.seqs[i];
        auto hi = hashnames.find(name);
        if (hi == nullptr) {
            throw std::invalid_argument(
                    strformat("Sequence %s from constraints file is not in the alignment", name.c_str()));
        }
        size_t iSeqNonunique = *hi;
        assert(iSeqNonunique >= 0 && iSeqNonunique < unique.alnToUniq.size());
        size_t iSeqUnique = unique.alnToUniq[iSeqNonunique];
        assert(iSeqUnique >= 0 && iSeqUnique < unique.uniqueSeq.size());
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