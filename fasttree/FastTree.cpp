

#include "FastTree.h"

using namespace fasttree;

FastTree::FastTree(const Options &options) : options(options) {}

void FastTree::run(std::istream &in, std::ostream &out, std::ostream &log) {
    options.codesString = options.nCodes == 20 ? Constants::codesStringAA : Constants::codesStringNT;

    if (options.nCodes == 4 && options.matrixPrefix.size() == 0)
        options.useMatrix = false;        /* no default nucleotide matrix */

    if (options.slow && options.tophitsMult > 0) {
        options.tophitsMult = 0.0;
    }

    if (!options.make_matrix) {        /* Report settings */
        std::string tophitString = "no";
        std::string tophitsCloseStr = "default";
/*





        char tophitString[100] = "no";
        char tophitsCloseStr[100] = "default";
        if (tophitsClose > 0) sprintf(tophitsCloseStr, "%.2f", tophitsClose);
        if (tophitsMult > 0)
            sprintf(tophitString, "%.2f*sqrtN close=%s refresh=%.2f",
                    tophitsMult, tophitsCloseStr, tophitsRefresh);
        char supportString[100] = "none";
        if (nBootstrap > 0) {
            if (MLnni != 0 || MLlen)
                sprintf(supportString, "SH-like %d", nBootstrap);
            else
                sprintf(supportString, "Local boot %d", nBootstrap);
        }
        char nniString[100] = "(no NNI)";
        if (nni > 0)
            sprintf(nniString, "+NNI (%d rounds)", nni);
        if (nni == -1)
            strcpy(nniString, "+NNI");
        char sprString[100] = "(no SPR)";
        if (spr > 0)
            sprintf(sprString, "+SPR (%d rounds range %d)", spr, maxSPRLength);
        char mlnniString[100] = "(no ML-NNI)";
        if (MLnni > 0)
            sprintf(mlnniString, "+ML-NNI (%d rounds)", MLnni);
        else if (MLnni == -1)
            sprintf(mlnniString, "+ML-NNI");
        else if (MLlen)
            sprintf(mlnniString, "+ML branch lengths");
        if ((MLlen || MLnni != 0) && !exactML)
            strcat(mlnniString, " approx");
        if (MLnni != 0)
            sprintf(mlnniString + strlen(mlnniString), " opt-each=%d", mlAccuracy);

        for (i = 0; i < nFPs; i++) {
            FILE *fp = fps[i];
            fprintf(fp, "FastTree Version %s %s%s\nAlignment: %s",
                    FT_VERSION, SSE_STRING, OpenMPString(), fileName != NULL ? fileName : "standard input");
            if (nAlign > 1)
                fprintf(fp, " (%d alignments)", nAlign);
            fprintf(fp, "\n%s distances: %s Joins: %s Support: %s\n",
                    nCodes == 20 ? "Amino acid" : "Nucleotide",
                    matrixPrefix ? matrixPrefix : (useMatrix ? "BLOSUM45"
                                                             : (nCodes == 4 && logdist ? "Jukes-Cantor"
                                                                                       : "%different")),
                    bionj ? "weighted" : "balanced",
                    supportString);
            if (intreeFile == NULL)
                fprintf(fp, "Search: %s%s %s %s %s\nTopHits: %s\n",
                        slow ? "Exhaustive (slow)" : (fastest ? "Fastest" : "Normal"),
                        useTopHits2nd ? "+2nd" : "",
                        nniString, sprString, mlnniString,
                        tophitString);
            else
                fprintf(fp, "Start at tree from %s %s %s\n", intreeFile, nniString, sprString);

            if (MLnni != 0 || MLlen) {
                fprintf(fp, "ML Model: %s,",
                        (nCodes == 4) ?
                        (bUseGtr ? "Generalized Time-Reversible" : "Jukes-Cantor") :
                        (bUseLg ? "Le-Gascuel 2008" : (bUseWag ? "Whelan-And-Goldman" : "Jones-Taylor-Thorton"))

                );
                if (nRateCats == 1)
                    fprintf(fp, " No rate variation across sites");
                else
                    fprintf(fp, " CAT approximation with %d rate categories", nRateCats);
                fprintf(fp, "\n");
                if (nCodes == 4 && bUseGtrRates)
                    fprintf(fp, "GTR rates(ac ag at cg ct gt) %.4f %.4f %.4f %.4f %.4f %.4f\n",
                            gtrrates[0], gtrrates[1], gtrrates[2], gtrrates[3], gtrrates[4], gtrrates[5]);
                if (nCodes == 4 && bUseGtrFreq)
                    fprintf(fp, "GTR frequencies(A C G T) %.4f %.4f %.4f %.4f\n",
                            gtrfreq[0], gtrfreq[1], gtrfreq[2], gtrfreq[3]);
            }
            if (constraintsFile != NULL)
                fprintf(fp, "Constraints: %s Weight: %.3f\n", constraintsFile, constraintWeight);
            if (pseudoWeight > 0)
                fprintf(fp, "Pseudocount weight for comparing sequences with little overlap: %.3lf\n", pseudoWeight);
            fflush(fp);
        }*/
    }


}
