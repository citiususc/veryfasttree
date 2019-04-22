
#include "FastTree.h"
#include "Utils.h"
#include "Debug.h"

using namespace fasttree;

FastTree::FastTree(const Options &options) : options(options) {}

void FastTree::prepare(std::istream &in, std::ostream &log) {
    options.codesString = options.nCodes == 20 ? Constants::codesStringAA : Constants::codesStringNT;

    if (options.nCodes == 4 && options.matrixPrefix.size() == 0)
        options.useMatrix = false;        /* no default nucleotide matrix */

    if (options.slow && options.tophitsMult > 0) {
        options.tophitsMult = 0.0;
    }

    if (!options.make_matrix) {        /* Report settings */
        std::string tophitString = "no";
        std::string tophitsCloseStr = "default";
        if (options.tophitsClose > 0) {
            tophitsCloseStr = strformat("%.2f", options.tophitsClose);
        }
        if (options.tophitsMult > 0) {
            tophitString = strformat("%.2f*sqrtN close=%s refresh=%.2f", options.tophitsMult, tophitsCloseStr,
                                     options.tophitsRefresh);
        }

        std::string supportString = "none";
        if (options.nBootstrap > 0) {
            if (options.MLnni != 0 || options.MLlen)
                supportString = strformat("SH-like %d", options.nBootstrap);
            else
                supportString = strformat("Local boot %d", options.nBootstrap);
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

        log << "FastTree Version " << Constants::version << " " << Constants::compileFlags << std::endl;
        log << "Alignment: " << (&in == &std::cin ? "standard input" : options.inFileName);

        if (options.nAlign > 1) {
            log << strformat(" (%d alignments)", options.nAlign) << std::endl;
        }
        log << (options.nCodes == 20 ? "Amino acid" : "Nucleotide ");
        log << "distances: " << (options.matrixPrefix.size() > 0 ? options.matrixPrefix :
                                 (options.useMatrix ? "BLOSUM45" : (options.nCodes == 4 && options.logdist
                                                                    ? "Jukes-Cantor"
                                                                    : "%different")));
        log << "Joins: " << (options.bionj ? "weighted" : "balanced");
        log << "Support: " << supportString << std::endl;


        if (options.intreeFile.size() == 0) {
            log << "Search: ";
            log << (options.slow ? "Exhaustive (slow)" : (options.fastest ? "Fastest" : "Normal"));
            log << (options.useTopHits2nd ? "+2nd" : "") << "";
            log << nniString << " " << sprString << " " << mlnniString;
            log << "TopHits: " << tophitString << std::endl;

        } else {
            log << "Start at tree from ";
            log << options.intreeFile << " ";
            log << nniString << " ";
            log << sprString << " ";
            log << std::endl;
        }

        if (options.MLnni != 0 || options.MLlen) {
            log << "ML Model: ",
                    log << ((options.nCodes == 4) ?
                            (options.bUseGtr ? "Generalized Time-Reversible" : "Jukes-Cantor") :
                            (options.bUseLg ? "Le-Gascuel 2008" : (options.bUseWag ? "Whelan-And-Goldman"
                                                                                   : "Jones-Taylor-Thorton")));
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
        if (options.constraintsFile.size() > 0) {
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

void FastTree::run(std::istream &in, std::ostream &out, std::ostream &log) {
    FastTreeImpl<double,DebugAllocator,DebugAllocator> impl(options);


}
