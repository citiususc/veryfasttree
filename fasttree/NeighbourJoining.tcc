
#ifndef FASTTREE_NEIGHBOURJOINING_TCC
#define FASTTREE_NEIGHBOURJOINING_TCC

#include "NeighbourJoining.h"
#include "assert.h"

#define AbsNeighbourJoining(...) \
template<typename Precision> \
__VA_ARGS__ fasttree::NeighbourJoining<Precision>


AbsNeighbourJoining()::Profile::Profile() {}

AbsNeighbourJoining()::Profile::Profile(size_t nPos, size_t nConstraints) {
    weights.resize(nPos);
    codes.resize(nPos);
    if (nConstraints == 0) {
        nOn.resize(nConstraints);
        nOff.resize(nConstraints);
    }
}

AbsNeighbourJoining()::Rates::Rates(size_t nRateCategories, size_t nPos) {
    assert(nRateCategories >= 0);
    if (nRateCategories > 0) {
        rates.resize(nRateCategories, 1.0);
        ratecat.resize(nPos, 0);
    }
}

AbsNeighbourJoining()::
NeighbourJoining(const Options &options, std::ostream &log, const std::vector<std::string> &seqs, size_t nPos,
                 const std::vector<std::string> &constraintSeqs,
                 DistanceMatrix <Precision> &distanceMatrix,
                 TransitionMatrix <Precision> &transmat) : log(log),
                                                           options(options),
                                                           seqs(seqs),
                                                           distanceMatrix(distanceMatrix),
                                                           transmat(transmat),
                                                           constraintSeqs(constraintSeqs),
                                                           outprofile(nPos, constraintSeqs.size()),
                                                           rates(1, nPos) {
    this->root = -1;
    this->maxnode = seqs.size();
    this->nPos = nPos;
    this->maxnodes = 2 * seqs.size();

    profiles.resize(maxnodes, Profile(nPos, constraintSeqs.size()));

    size_t counts[256] = {}; /*Array of zeros*/

    for (int c = 0; c < 256; c++) {
        charToCode[c] = static_cast<uint8_t>(options.nCodes);
    }
    for (int i = 0; options.codesString[i]; i++) {
        charToCode[static_cast<int>(options.codesString[i])] = static_cast<uint8_t>(i);
    }
    charToCode['-'] = NOCODE;

    for (size_t i = 0; i < seqs.size(); i++) {//TODO OpenMP
        seqToProfile(i, profiles[i], seqs[i], constraintSeqs[i]);
    }

    size_t totCount = 0;
    for (int i = 0; i < 256; i++)
        totCount += counts[i];

    /* warnings about unknown characters */
    for (int i = 0; i < 256; i++) {
        if (counts[i] == 0 || i == '-') {
            continue;
        }

        if (options.codesString.find(static_cast<uint8_t>(i)) == std::string::npos) {
            log << strformat("Ignored unknown character %c (seen %lu times)", i, counts[i]) << std::endl;
        }
    }

    /* warnings about the counts */
    double fACGTUN = (counts['A'] + counts['C'] + counts['G'] + counts['T'] + counts['U'] + counts['N'])
                     / (double) (totCount - counts['-'] - counts['.']);
    if (options.nCodes == 4 && fACGTUN < 0.9) {
        log << strformat("WARNING! ONLY %.1f%% NUCLEOTIDE CHARACTERS -- IS THIS REALLY A NUCLEOTIDE ALIGNMENT?",
                         100.0 * fACGTUN) << std::endl;
    } else if (options.nCodes == 20 && fACGTUN >= 0.9) {
        log << strformat("WARNING! %.1f%% NUCLEOTIDE CHARACTERS -- IS THIS REALLY A PROTEIN ALIGNMENT?",
                         100.0 * fACGTUN) << std::endl;
    }

    if (options.verbose > 10) {
        log << "Made sequence profiles" << std::endl;
    }

    /* profiles from nSeq to maxnodes not yet exists */

    outProfile(outprofile);

    if (options.verbose > 10) {
        log << "Made out-profile" << std::endl;
    }

    totdiam = 0.0;

    diameter.resize(maxnodes, 0);
    varDiameter.resize(maxnodes, 0);
    selfdist.resize(maxnodes, 0);
    selfweight.resize(maxnodes);
    for (size_t i = 0; i < seqs.size(); i++) {//TODO OpenMP
        selfweight[i] = nPos - profiles[i].nGaps;
        assert(selfweight[i] == (nPos - nGaps(i)));
    }

    outDistances.resize(maxnodes);
    nOutDistActive.resize(maxnodes, seqs.size() * 10); /* unreasonably high value */
    // parent.empty()        /* so SetOutDistance ignores it */
    for (size_t i = 0; i < seqs.size(); i++) {
        setOutDistance(i, /*nActive*/seqs.size());
    }

    if (options.verbose > 2) {
        for (size_t i = 0; i < 4 && i < seqs.size(); i++) {
            log << strformat("Node %d outdist %f", i, outDistances[i]) << std::endl;
        }
    }

    parent.resize(maxnodes, -1);
    branchlength.resize(maxnodes, 0); /* distance to parent */
    support.resize(maxnodes, -1.0);
    child.resize(maxnodes);
}

AbsNeighbourJoining(void)::
seqToProfile(size_t i, Profile &profile, const std::string &seq, const std::string &constraintSeq) {
//TODO

}

AbsNeighbourJoining(void)::outProfile(Profile &out) {
//TODO
}


AbsNeighbourJoining(void)::setOutDistance(size_t i, size_t nActive) {
//TODO
}

AbsNeighbourJoining(size_t)::nGaps(size_t i) {
    assert(i < seqs.size());
    int nGaps = 0;
    for (size_t p = 0; p < nPos; p++) {
        if (profiles[i].codes[p] == NOCODE) {
            nGaps++;
        }
    }
    return nGaps;
}

#endif