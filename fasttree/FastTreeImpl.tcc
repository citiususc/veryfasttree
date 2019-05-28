
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

    if (options.constraintsFile.empty()) {
        fpConstraints.setstate(std::ios_base::badbit);
    } else {
        fpConstraints.open(options.constraintsFile);
        if (fpConstraints.fail()) {
            throw std::invalid_argument("Cannot read " + options.constraintsFile);
        }
    }

    if (options.intreeFile.empty()) {
        fpInTree.setstate(std::ios_base::badbit);
    } else {
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
            log << strformat("Read %zd sequences, %ld positions", aln.seqs.size(), aln.nPos) << std::endl;
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
                    log << strformat("alignment #%ld in ", iAln + 1) << options.constraintsFile << std::endl;
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
            int64_t nSeq = (int64_t) aln.seqs.size();
            int64_t unSeq = (int64_t) unique.uniqueSeq.size();
            int64_t unConstraints = (int64_t) uniqConstraints.size();

            if (options.verbose > 2) {
                log << strformat("read %s seqs %ld (%ld unique) positions %ld nameLast %s seqLast %s",
                                 options.inFileName.empty() ? "standard input" : options.inFileName.c_str(),
                                 nSeq, unSeq, aln.nPos,
                                 aln.names[aln.seqs.size() - 1].c_str(),
                                 aln.seqs[aln.seqs.size() - 1].c_str()) << std::endl;
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

            int64_t nniToDo = options.nni == -1 ? (int64_t) (0.5 + 4.0 * std::log(unSeq) / std::log(2)) :
                              options.nni;
            int64_t sprRemaining = options.spr;
            int64_t MLnniToDo = (options.MLnni != -1) ? options.MLnni :
                                (int64_t) (0.5 + 2.0 * std::log(unSeq) / std::log(2));
            if (options.verbose > 0) {
                if (!fpInTree) {
                    log << strformat("Initial topology in %.2f seconds", progressReport.clockDiff()) << std::endl;
                }
                if (options.spr > 0 || nniToDo > 0 || MLnniToDo > 0) {
                    log << strformat("Refining topology: %ld rounds ME-NNIs, %d rounds ME-SPRs, %ld rounds ML-NNIs",
                                     nniToDo, options.spr, MLnniToDo) << std::endl;
                }
            }

            if (nniToDo > 0) {
                bool bConverged = false;
                std::vector<typename decltype(nj)::NNIStats> nni_stats;
                nj.initNNIStats(nni_stats);

                for (int64_t i = 0; i < nniToDo; i++) {
                    double maxDelta;
                    if (!bConverged) {
                        int64_t nChange = nj.DoNNI(i, nniToDo,/*use ml*/false, nni_stats, maxDelta);

                        nj.logTree("ME_NNI%ld", i + 1, aln.names, unique);
                        if (nChange == 0) {
                            bConverged = true;
                            if (options.verbose > 1) {
                                log << strformat("Min_evolution NNIs converged at round %ld -- skipping some rounds",
                                                 i + 1) << std::endl;
                            }
                            if (!options.logFileName.empty()) {
                                log << strformat("Min_evolution NNIs converged at round %ld -- skipping some rounds",
                                                 i + 1) << std::endl;
                            }
                        }
                    }

                    /* Interleave SPRs with NNIs (typically 1/3rd NNI, SPR, 1/3rd NNI, SPR, 1/3rd NNI */
                    if (sprRemaining > 0 &&
                        (nniToDo / (options.spr + 1) > 0 &&
                         ((i + 1) % (nniToDo / (options.spr + 1))) == 0)) {
                        nj.SPR(options.maxSPRLength, options.spr - sprRemaining, options.spr);
                        nj.logTree("ME_SPR%ld", options.spr - sprRemaining + 1, aln.names, unique);
                        sprRemaining--;
                        /* Restart the NNIs -- set all ages to 0, etc. */
                        bConverged = false;
                        nj.initNNIStats(nni_stats);
                    }
                }
            }
            while (sprRemaining > 0) {    /* do any remaining SPR rounds */
                nj.SPR(options.maxSPRLength, options.spr - sprRemaining, options.spr);
                nj.logTree("ME_SPR%ld", options.spr - sprRemaining + 1, aln.names, unique);
                sprRemaining--;
            }

            /* In minimum-evolution mode, update branch lengths, even if no NNIs or SPRs,
           so that they are log-corrected, do not include penalties from constraints,
           and avoid errors due to approximation of out-distances.
           If doing maximum-likelihood NNIs, then we'll also use these
           to get estimates of starting distances for quartets, etc.
          */
            nj.updateBranchLengths();
            nj.logTree("ME_Lengths", 0, aln.names, unique);

            double total_len = nj.totalLen();

            if (options.verbose > 0) {
                log << strformat("Total branch-length %.3f after %.2f sec",
                                 total_len, progressReport.clockDiff()) << std::endl;
            }

            typename decltype(nj)::SplitCount splitcount = {0, 0, 0, 0, 0.0, 0.0};

            if (MLnniToDo > 0 || options.MLlen) {
                bool warn_len =
                        total_len / nj.getMaxnode() < 0.001 && options.MLMinBranchLengthTolerance > 1.0 / aln.nPos;
                bool warn = warn_len || (total_len / nj.getMaxnode() < 0.001 && aln.nPos >= 10000);

                if (warn) {
                    log << std::endl;
                    log << "WARNING! This alignment consists of closely-related and very-long sequences." << std::endl;
                }
                if (warn_len) {
                    log << "This version of FastTree may not report reasonable branch lengths!" << std::endl;
                    if (options.doublePrecision) {
                        log << "Consider changing MLMinBranchLengthTolerance." << std::endl;
                    } else {
                        log << "Consider use FastTree with -double-precision." << std::endl;
                    }
                    log << "For more information, visit" << std::endl;
                    log << "http://www.microbesonline.org/fasttree/#BranchLen" << std::endl;

                }
                if (warn) {
                    log << "WARNING! FastTree (or other standard maximum-likelihood tools)" << std::endl;
                    log << "may not be appropriate for aligments of very closely-related sequences" << std::endl;
                    log << "like this one, as FastTree does not account for recombination or gene conversion"
                        << std::endl << std::endl;

                }

                /* Do maximum-likelihood computations */
                DistanceMatrix <Precision, op_t::ALIGNMENT> tmatAsDist;
                /* Convert profiles to use the transition matrix */
                transMatToDistanceMat(tmatAsDist);
                nj.recomputeProfiles(tmatAsDist);

                double lastloglk = -1e20;
                std::vector<typename decltype(nj)::NNIStats> nni_stats;
                nj.initNNIStats(nni_stats);
                bool resetGtr = options.nCodes == 4 && options.bUseGtr && !options.bUseGtrRates;

                if (options.MLlen) {
                    int64_t maxRound = (int64_t) (0.5 + std::log(unSeq) / std::log(2));
                    double dLastLogLk = -1e20;
                    for (int64_t iRound = 1; iRound <= maxRound; iRound++) {
                        auto &branchlength = nj.getBranchlength();
                        std::vector<numeric_t, typename op_t::Allocator> oldlength(branchlength.begin(),
                                                                                   branchlength.begin() +
                                                                                   nj.getMaxnode());

                        nj.optimizeAllBranchLengths();
                        nj.logTree("ML_Lengths", iRound, aln.names, unique);
                        double dMaxChange = 0; /* biggest change in branch length */
                        for (int64_t node = 0; node < nj.getMaxnode(); node++) {
                            double d = std::fabs(oldlength[node] - branchlength[node]);
                            if (dMaxChange < d) {
                                dMaxChange = d;
                            }
                        }

                        double loglk = nj.treeLogLk(/*site_likelihoods*/nullptr);
                        bool bConverged = iRound > 1 &&
                                          (dMaxChange < 0.001 || loglk < (dLastLogLk + Constants::treeLogLkDelta));
                        if (options.verbose) {
                            log << strformat("%ld rounds ML lengths: LogLk %s= %.3lf Max-change %.4lf%s Time %.2f",
                                             iRound,
                                             options.exactML || options.nCodes != 20 ? "" : "~",
                                             loglk,
                                             dMaxChange,
                                             bConverged ? " (converged)" : "",
                                             progressReport.clockDiff()) << std::endl;
                        }
                        if (!options.logFileName.empty()) {
                            log << strformat("TreeLogLk\tLength%ld\t%.4lf\tMaxChange\t%.4lf",
                                             iRound, loglk, dMaxChange) << std::endl;
                        }
                        if (iRound == 1) {
                            if (resetGtr) {
                                nj.setMLGtr(options.bUseGtrFreq ? options.gtrfreq : nullptr);
                            }
                            nj.setMLRates();
                            nj.logMLRates();
                        }
                        if (bConverged) {
                            break;
                        }
                    }
                }

                if (MLnniToDo > 0) {
                    /* This may help us converge faster, and is fast */
                    nj.optimizeAllBranchLengths();
                    nj.logTree("ML_Lengths%ld", 1, aln.names, unique);
                }

                double maxDelta;
                bool bConverged = false;
                for (int64_t iMLnni = 0; iMLnni < MLnniToDo; iMLnni++) {
                    int64_t changes = nj.DoNNI(iMLnni, MLnniToDo, /*use ml*/true, nni_stats, maxDelta);
                    nj.logTree("ML_NNI%ld", iMLnni + 1, aln.names, unique);
                    double loglk = nj.treeLogLk(/*site_likelihoods*/nullptr);
                    bool bConvergedHere = (iMLnni > 0) &&
                                          ((loglk < lastloglk + Constants::treeLogLkDelta) ||
                                           maxDelta < Constants::treeLogLkDelta);
                    if (options.verbose) {
                        log << strformat("ML-NNI round %ld: LogLk %s= %.3f NNIs %ld max delta %.2f Time %.2f%s",
                                         iMLnni + 1,
                                         options.exactML || options.nCodes != 20 ? "" : "~",
                                         loglk, changes, maxDelta, progressReport.clockDiff(),
                                         bConverged ? " (final)" : "") << std::endl;
                    }
                    if (!options.logFileName.empty())
                        log << strformat("TreeLogLk\tML_NNI%ld\t%.4lf\tMaxChange\t%.4lf",
                                         iMLnni + 1, loglk, maxDelta) << std::endl;
                    if (bConverged) {
                        break;        /* we did our extra round */
                    }
                    if (bConvergedHere) {
                        bConverged = true;
                    }
                    if (bConverged || iMLnni == MLnniToDo - 2) {
                        /* last round uses high-accuracy seettings -- reset NNI stats to tone down heuristics */
                        nni_stats.clear();
                        nj.initNNIStats(nni_stats);
                        if (options.verbose) {
                            log << "Turning off heuristics for final round of ML NNIs"
                                << (bConvergedHere ? " (converged)" : "") << std::endl;
                        }

                    }
                    lastloglk = loglk;
                    if (iMLnni == 0 && nj.getRateCategories() == 1) {
                        if (resetGtr) {
                            nj.setMLGtr(options.bUseGtrFreq ? options.gtrfreq : nullptr);
                        }
                        nj.setMLRates();
                        nj.logMLRates();
                    }
                }
                nni_stats.clear();
                nni_stats.reserve(0);


                /* This does not take long and improves the results */
                if (MLnniToDo > 0) {
                    nj.optimizeAllBranchLengths();
                    nj.logTree("ML_Lengths%ld", 2, aln.names, unique);
                    if (options.verbose || !options.logFileName.empty()) {
                        double loglk = nj.treeLogLk(/*site_likelihoods*/nullptr);
                        if (options.verbose)
                            log << strformat("Optimize all lengths: LogLk %s= %.3f Time %.2f",
                                             options.exactML || options.nCodes != 20 ? "" : "~",
                                             loglk,
                                             progressReport.clockDiff()) << std::endl;
                        if (!options.logFileName.empty()) {
                            log << strformat("TreeLogLk\tML_Lengths%d\t%.4f", 2, loglk) << std::endl;
                        }
                    }
                }

                /* Count bad splits and compute SH-like supports if desired */
                if ((MLnniToDo > 0 && !options.fastest) || options.nBootstrap > 0)
                    nj.testSplitsML(splitcount);

                /* Compute gamma-based likelihood? */
                if (options.gammaLogLk && options.nRateCats > 1) {
                    nj.branchlengthScale();
                }
            } else {
                /* Minimum evolution supports */
                nj.testSplitsMinEvo(splitcount);
                if (options.nBootstrap > 0) {
                    nj.reliabilityNJ();
                }
            }

            log << strformat("Total time: %.2f seconds Unique: %ld/%ld Bad splits: %ld/%ld",
                             progressReport.clockDiff(),
                             unSeq, nSeq,
                             splitcount.nBadSplits, splitcount.nSplits);
            if (splitcount.dWorstDeltaUnconstrained > 0) {
                log << strformat(" Worst %sdelta-%s %.3f",
                                 !uniqConstraints.empty() ? "unconstrained " : "",
                                 (MLnniToDo > 0 || options.MLlen) ? "LogLk" : "Len",
                                 splitcount.dWorstDeltaUnconstrained);
            }
            log << std::endl;
            if (unSeq > 3 && unConstraints > 0) {
                log << strformat("Violating constraints: %ld both bad: %ld",
                                 splitcount.nConstraintViolations, splitcount.nBadBoth);
                if (splitcount.dWorstDeltaConstrained > 0) {
                    log << strformat(" Worst delta-%s due to constraints: %.3f",
                                     (MLnniToDo > 0 || options.MLlen) ? "LogLk" : "Len",
                                     splitcount.dWorstDeltaConstrained);
                }
                log << std::endl;
            }
            if (options.threads > 1) {

            } else if (options.verbose > 1 || !options.logFileName.empty()) {
                double dN2 = unSeq * (double) unSeq;
                log << strformat("Dist/N**2: by-profile %.3f (out %.3f) by-leaf %.3f avg-prof %.3f",
                                 options.debug.profileOps / dN2,
                                 options.debug.outprofileOps / dN2,
                                 options.debug.seqOps / dN2,
                                 options.debug.profileAvgOps / dN2) << std::endl;

                if (options.debug.nCloseUsed > 0 || options.debug.nClose2Used > 0 ||
                    options.debug.nRefreshTopHits > 0) {
                    log << strformat("Top hits: close neighbors %ld/%ld 2nd-level %ld refreshes %ld",
                                     options.debug.nCloseUsed, unSeq, options.debug.nClose2Used,
                                     options.debug.nRefreshTopHits);
                }
                if (!options.slow) {
                    log << strformat(" Hill-climb: %ld Update-best: %ld", options.debug.nHillBetter,
                                     options.debug.nVisibleUpdate) << std::endl;
                }
                if (nniToDo > 0 || options.spr > 0 || MLnniToDo > 0) {
                    log << strformat("NNI: %ld SPR: %ld ML-NNI: %ld", options.debug.nNNI, options.debug.nSPR,
                                     options.debug.nML_NNI) << std::endl;
                }
                if (MLnniToDo > 0) {
                    log << strformat("Max-lk operations: lk %ld posterior %ld", options.debug.nLkCompute,
                                     options.debug.nPosteriorCompute);
                    if (options.debug.nAAPosteriorExact > 0 || options.debug.nAAPosteriorRough > 0) {
                        log << strformat(" approximate-posteriors %.2f%%",
                                         (100.0 * options.debug.nAAPosteriorRough) /
                                         (double) (options.debug.nAAPosteriorExact + options.debug.nAAPosteriorRough));
                    }
                    if (options.mlAccuracy < 2) {
                        log << strformat(" star-only %ld", options.debug.nStarTests);
                    }
                    log << std::endl;
                }
            }


            nj.printNJ(output, aln.names, unique, /*support*/options.nBootstrap > 0);
            log << "TreeCompleted" << std::endl;
        }/* end build tree */
    }/* end loop over alignments */
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

AbsFastTreeImpl(void)::transMatToDistanceMat(DistanceMatrix <Precision, op_t::ALIGNMENT> &dmat) {
    if (!transmat) {
        return;
    }

    for (int64_t i = 0; i < options.nCodes; i++) {
        for (int64_t j = 0; j < options.nCodes; j++) {
            dmat.distances[i][j] = 0;    /* never actually used */
            dmat.eigeninv[i][j] = transmat.eigeninv[i][j];
            dmat.codeFreq[i][j] = transmat.codeFreq[i][j];
        }
    }
    /* eigentot . rotated-vector is the total frequency of the unrotated vector
       (used to normalize in NormalizeFreq()
       For transition matrices, we rotate by transpose of eigenvectors, so
       we need to multiply by the inverse matrix by 1....1 to get this vector,
       or in other words, sum the columns
    */
    for (int64_t i = 0; i < options.nCodes; i++) {
        dmat.eigentot[i] = 0.0;
        for (int64_t j = 0; j < options.nCodes; j++)
            dmat.eigentot[i] += transmat.eigeninv[i][j];
    }
    dmat.setted = true;
}


#endif