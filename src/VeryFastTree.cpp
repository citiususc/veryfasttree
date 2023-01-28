
#include "VeryFastTree.h"
#include "Utils.h"
#include "Debug.h"
#include <omp.h>
#include "operations/BasicOperations.h"

#if (defined __SSE2__) || (defined __AVX__)

#include "operations/SSE128Operations.h"

#endif

#ifdef __AVX__

#include "operations/AVX256Operations.h"

#endif

#ifdef __AVX512F__

#include "operations/AVX512Operations.h"

#endif

using namespace veryfasttree;

VeryFastTree::VeryFastTree(const Options &options) : options(options) {}

void VeryFastTree::settings(std::ostream &log) {
    if (options.extension == "AUTO") {
        options.extension = "NONE";
        #ifdef __SSE2__
        options.extension = "SSE3";
        #endif
        #ifdef __AVX__
        if (options.doublePrecision) {
            options.extension = "AVX2";
        }
        #endif
    }

    if (options.nCodes == 4){
        if (options.extension == "AVX512") {
            options.extension = "AVX2";
        }
        if(!options.doublePrecision && options.extension == "AVX2"){
            options.extension = "SSE3";
        }
    }

    if (options.threads > 1){
        options.verbose = 1;
    }

    if(options.useTopHits2nd && options.threads > 1){
        options.useTopHits2nd = false;
        log << "Warning: 2nd-level will be ignored because using at least 2 threads results in the best performance "
               "without marginal reductions in tree quality";
    }


    options.codesString = options.nCodes == 20 ? Constants::codesStringAA : Constants::codesStringNT;

    if (options.nCodes == 4 && options.matrixPrefix.empty()) {
        options.useMatrix = false;        /* no default nucleotide matrix */
    }

    if (!options.transitionFile.empty() && options.nCodes != 20) {
        throw std::invalid_argument("The -trans option is only supported for amino acid alignments");
    }

    if (!options.transitionFile.empty()) {
        log << "Warning: custom matrices may create numerical problems for single-precision VeryFastTree.\n"
               "You may want to use with -double-precision";
    }

    if (options.slow && options.fastest) {
        throw std::invalid_argument("Cannot be both slow and fastest");
    }

    if (options.slow && options.tophitsMult > 0) {
        options.tophitsMult = 0.0;
    }

    if (options.doublePrecision) {
        options.MLMinBranchLengthTolerance = Constants::MLMinBranchLengthToleranceDouble;
        options.MLFTolBranchLength = Constants::MLFTolBranchLengthDouble;
        options.MLMinBranchLength = Constants::MLMinBranchLengthDouble;
        options.MLMinRelBranchLength = Constants::MLMinRelBranchLengthDouble;
        options.fPostTotalTolerance = Constants::fPostTotalToleranceDouble;
    } else {
        options.MLMinBranchLengthTolerance = Constants::MLMinBranchLengthToleranceFloat;
        options.MLFTolBranchLength = Constants::MLFTolBranchLengthFloat;
        options.MLMinBranchLength = Constants::MLMinBranchLengthFloat;
        options.MLMinRelBranchLength = Constants::MLMinRelBranchLengthFloat;
        options.fPostTotalTolerance = Constants::fPostTotalToleranceFloat;
    }

    if (!options.makeMatrix) {        /* Report settings */
        std::string tophitString = "no";
        std::string tophitsCloseStr = "default";
        if (options.tophitsClose > 0) {
            tophitsCloseStr = strformat("%.2f", options.tophitsClose);
        }
        if (options.tophitsMult > 0) {
            tophitString = strformat("%.2f*sqrtN close=%s refresh=%.2f", options.tophitsMult, tophitsCloseStr.c_str(),
                                     options.tophitsRefresh);
        }

        std::string supportString = "none";
        if (options.nBootstrap > 0) {
            if (options.MLnni != 0 || options.MLlen) {
                supportString = strformat("SH-like %d", options.nBootstrap);
            } else {
                supportString = strformat("Local boot %d", options.nBootstrap);
            }
        }
        std::string nniString = "(no NNI)";
        if (options.nni > 0) {
            nniString = strformat("+NNI (%d rounds)", options.nni);
        }
        if (options.nni == -1) {
            nniString = "+NNI";
        }
        std::string sprString = "(no SPR)";
        if (options.spr > 0) {
            sprString = strformat("+SPR (%d rounds range %d)", options.spr, options.maxSPRLength);
        }
        std::string mlnniString = "(no ML-NNI)";
        if (options.MLnni > 0) {
            mlnniString = strformat("+ML-NNI (%d rounds)", options.MLnni);
        } else if (options.MLnni == -1) {
            mlnniString = "+ML-NNI";
        } else if (options.MLlen) {
            mlnniString = "+ML branch lengths";
        }
        if ((options.MLlen || options.MLnni != 0) && !options.exactML) {
            mlnniString += " approx";
        }
        if (options.MLnni != 0) {
            mlnniString += strformat(" opt-each=%d", options.mlAccuracy);
        }

        log << "VeryFastTree Version " << Constants::version << " " << Constants::compileFlags;
        if (options.extension != "NONE") {
            log << " with " << options.extension;
        }
        if (options.threads > 1) {
            log << " using " << "threads(" << options.threads << ") ";
            log << "level " << options.threadsLevel;
            if (options.deterministic) {
                log << " deterministic";
            } else {
                log << " NO deterministic";
            }
        }
        log << std::endl;
        log << "Alignment: " << (options.inFileName.empty() ? "standard input" : options.inFileName);

        if (options.nAlign > 1) {
            log << strformat(" (%ld alignments)", options.nAlign) << std::endl;
        }
        log << std::endl << strformat(
                "%s distances: %s Joins: %s Support: %s",
                options.nCodes == 20 ? "Amino acid" : "Nucleotide",
                !options.matrixPrefix.empty() ? options.matrixPrefix.c_str() :
                (options.useMatrix ? "BLOSUM45" : (options.nCodes == 4 && options.logdist
                                                   ? "Jukes-Cantor"
                                                   : "%different")),
                options.bionj ? "weighted" : "balanced",
                supportString.c_str()
        ) << std::endl;

        if (options.intreeFile.empty()) {
            log << strformat(
                    "Search: %s%s %s %s %s",
                    (options.slow ? "Exhaustive (slow)" : (options.fastest ? "Fastest" : "Normal")),
                    (options.useTopHits2nd ? "+2nd" : ""),
                    nniString.c_str(), sprString.c_str(), mlnniString.c_str()
            ) << std::endl;
            log << strformat("TopHits: %s", tophitString.c_str()) << std::endl;
        } else {
            log << strformat("Start at tree from %s %s %s",
                             options.intreeFile.c_str(), nniString.c_str(), sprString.c_str()

            ) << std::endl;
        }

        if (options.MLnni != 0 || options.MLlen) {
            log << "ML Model: ",
                    log << ((options.nCodes == 4) ?
                            (options.bUseGtr ? "Generalized Time-Reversible" : "Jukes-Cantor") :
                            (!options.transitionFile.empty() ? options.transitionFile :
                             (options.bUseLg ? "Le-Gascuel 2008" : (options.bUseWag ? "Whelan-And-Goldman"
                                                                                    : "Jones-Taylor-Thorton"))));
            log << ",";
            if (options.nRateCats == 1) {
                log << " No rate variation across sites";
            } else {
                log << strformat(" CAT approximation with %d rate categories", options.nRateCats);
            }
            log << std::endl;
            if (options.nCodes == 4 && options.bUseGtrRates) {
                log << strformat("GTR rates(ac ag at cg ct gt) %.4f %.4f %.4f %.4f %.4f %.4f",
                                 options.gtrrates[0], options.gtrrates[1], options.gtrrates[2],
                                 options.gtrrates[3], options.gtrrates[4], options.gtrrates[5]) << std::endl;
            }
            if (options.nCodes == 4 && options.bUseGtrFreq) {
                log << strformat("GTR frequencies(A C G T) %.4f %.4f %.4f %.4f",
                                 options.gtrfreq[0], options.gtrfreq[1], options.gtrfreq[2], options.gtrfreq[3])
                    << std::endl;
            }
        }
        if (!options.constraintsFile.empty()) {
            log << "Constraints: " << options.constraintsFile;
            log << strformat(" Weight: %.3f", options.constraintWeight) <<
                std::endl;
        }
        if (options.pseudoWeight > 0) {
            log << strformat("Pseudocount weight for comparing sequences with little overlap: %.3lf",
                             options.pseudoWeight) << std::endl;
        }
        log.flush();
    }
}

void VeryFastTree::configOpenMP() {
    omp_set_num_threads(options.threads);
}

void VeryFastTree::run(std::istream &in, std::ostream &out, std::ostream &log) {
    settings(log);
    configOpenMP();

    if (options.extension == "NONE") {
        if (options.doublePrecision) {
            VeyFastTreeImpl<double, BasicOperations>(options, in, out, log).run();
        } else {
            VeyFastTreeImpl<float, BasicOperations>(options, in, out, log).run();
        }
    }
            #if (defined __SSE2__) || (defined __AVX__)
    else if (options.extension == "SSE" || options.extension == "SSE3") {
        if (options.doublePrecision) {
            VeyFastTreeImpl<double, SSE128Operations>(options, in, out, log).run();
        } else {
            VeyFastTreeImpl<float, SSE128Operations>(options, in, out, log).run();
        }
    }
            #endif
            #ifdef __AVX__
    else if (options.extension == "AVX" || options.extension == "AVX2") {
        if (options.doublePrecision) {
            VeyFastTreeImpl<double, AVX256Operations>(options, in, out, log).run();
        } else {
            VeyFastTreeImpl<float, AVX256Operations>(options, in, out, log).run();
        }
    }
            #endif
            #ifdef __AVX512F__
    else if (options.extension == "AVX512") {
        if (options.doublePrecision) {
            VeyFastTreeImpl<double, AVX512Operations>(options, in, out, log).run();
        } else {
            VeyFastTreeImpl<float, AVX512Operations>(options, in, out, log).run();
        }
    }
            #endif
    else {
        throw std::invalid_argument("This version has not been compiled with " + options.extension + " extension");
    }

}
