
#ifndef FASTTREE_NEIGHBOURJOINING_TCC
#define FASTTREE_NEIGHBOURJOINING_TCC

#include "NeighbourJoining.h"
#include "assert.h"

#define AbsNeighbourJoining(...) \
template<typename Precision, template<typename> typename Operations> \
__VA_ARGS__ fasttree::NeighbourJoining<Precision,Operations>


AbsNeighbourJoining()::Profile::Profile() {}

AbsNeighbourJoining()::Profile::Profile(size_t nPos, size_t nConstraints) {
    weights.resize(nPos);
    codes.resize(nPos);
    if (nConstraints == 0) {
        nOn.resize(nConstraints, 0);
        nOff.resize(nConstraints, 0);
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
NeighbourJoining(Options &options, std::ostream &log, std::vector<std::string> &seqs, size_t nPos,
                 std::vector<std::string> &constraintSeqs,
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


    seqsToProfiles();
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
    for (size_t i = 0; i < seqs.size(); i++) {
        selfweight[i] = nPos - profiles[i].nGaps;
        assert(selfweight[i] == (nPos - nGaps(i)));
    }

    outDistances.resize(maxnodes);
    nOutDistActive.resize(maxnodes, seqs.size() * 10); /* unreasonably high value */
    // parent.empty()        /* so SetOutDistance ignores it */
    #pragma omp parallel for schedule(static)
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

AbsNeighbourJoining(void)::printDistances(std::vector<std::string> &names, std::ostream &out) {
    for (size_t i = 0; i < seqs.size(); i++) {
        std::cout << names[i];
        for (size_t j = 0; j < seqs.size(); j++) {
            Besthit hit;
            seqDist(profiles[i].codes, profiles[j].codes, hit);
            if (options.logdist) {
                hit.dist = logCorrect(hit.dist);
            }
            /* Make sure -0 prints as 0 */
            std::cout << strformat(" %f", hit.dist <= 0.0 ? 0.0 : hit.dist);
        }
        std::cout << std::endl;
    }
}

AbsNeighbourJoining(double)::logCorrect(double dist) {
    const double maxscore = 3.0;
    if (options.nCodes == 4 && !options.useMatrix) { /* Jukes-Cantor */
        dist = dist < 0.74 ? -0.75 * std::log(1.0 - dist * 4.0 / 3.0) : maxscore;
    } else {            /* scoredist-like */
        dist = dist < 0.99 ? -1.3 * std::log(1.0 - dist) : maxscore;
    }
    return (dist < maxscore ? dist : maxscore);
}

AbsNeighbourJoining(void)::seqsToProfiles() {
    profiles.resize(maxnodes, Profile(nPos, constraintSeqs.size()));
    uint8_t charToCode[256];
    size_t counts[256] = {}; /*Array of zeros*/

    for (int c = 0; c < 256; c++) {
        charToCode[c] = static_cast<uint8_t>(options.nCodes);
    }
    for (int i = 0; options.codesString[i]; i++) {
        charToCode[static_cast<int>(options.codesString[i])] = static_cast<uint8_t>(i);
    }
    charToCode['-'] = NOCODE;

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < seqs.size(); i++) {
        auto &seq = seqs[i];
        auto &profile = profiles[i];
        for (size_t j = 0; j < nPos; j++) {
            auto character = (uint8_t) seq[j];
            counts[character]++;
            auto c = charToCode[character];
            if (options.verbose > 10 && j < 2) {
                log << strformat("pos %d char %c code %d", j, seq[j], c) << std::endl;
            }
            /* treat unknowns as gaps */
            if (c == options.nCodes || c == NOCODE) {
                profile.codes[j] = NOCODE;
                profile.weights[j] = 0.0;
                profile.nGaps++;
            } else {
                profile.codes[j] = c;
                profile.weights[j] = 1.0;
            }
        }
        if (!constraintSeqs.empty()) {
            auto& constraintSeq = constraintSeqs[i];
            for (size_t j = 0; j < constraintSeq.size(); j++) {
                if (constraintSeq[j] == '1') {
                    profile.nOn[j] = 1;
                } else if (constraintSeq[j] == '0') {
                    profile.nOff[j] = 1;
                } else if (constraintSeq[j] != '-') {
                    #pragma omp critical
                    {
                        log << strformat("Constraint characters in unique sequence %d replaced with gap: %c%d",
                                         i + 1, constraintSeq[j], j + 1) << std::endl;

                    };
                    /* For the benefit of ConstraintSequencePenalty -- this is a bit of a hack, as
                       this modifies the value read from the alignment
                    */
                    constraintSeq[j] = '-';
                }
            }
        }
    }

    size_t totCount = 0;
    for (int i = 0; i < 256; i++) {
        totCount += counts[i];
    }

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
}

AbsNeighbourJoining(void)::outProfile(Profile &out) {
    double inweight = 1.0 / (double) profiles.size();   /* The maximal output weight is 1.0 */

    /* First, set weights -- code is always NOCODE, prevent weight=0 */
    size_t nVectors = 0;
    for (size_t i = 0; i < nPos; i++) {
        out.weights[i] = 0;
        for (size_t in = 0; in < profiles.size(); in++) {
            out.weights[i] += profiles[in].weights[i] * inweight;
        }
        if (out.weights[i] <= 0) {
            out.weights[i] = 1e-20;
        } /* always store a vector */
        nVectors++;
        out.codes[i] = NOCODE;        /* outprofile is normally complicated */
    }

    /* Initialize the frequencies to 0 */
    out.vectors.reserve(nVectors * options.nCodes);
    out.vectors.resize(nVectors, 0);

    /* Add up the weights, going through each sequence in turn */
    #pragma omp parallel for schedule(static)
    for (size_t in = 0; in < profiles.size(); in++) {
        size_t iFreqOut = 0;
        size_t iFreqIn = 0;
        for (size_t i = 0; i < nPos; i++) {
            numeric_t *fIn = getFreq(profiles[in], i, iFreqIn);
            numeric_t *fOut = getFreq(out, i, iFreqOut);
            if (profiles[in].weights[i] > 0) {
                addToFreq(fOut, profiles[in].weights[i], profiles[in].codes[i], fIn);
            }
        }
        assert(iFreqOut == out.vectors.size());
        assert(iFreqIn == profiles[in].vectors.size());
    }

    /* And normalize the frequencies to sum to 1 */
    size_t iFreqOut = 0;
    for (size_t i = 0; i < nPos; i++) {
        numeric_t *fOut = getFreq(out, i, iFreqOut);
        if (fOut != nullptr) {
            normalizeFreq(fOut);
        }
    }
    assert(iFreqOut == out.vectors.size());
    if (options.verbose > 10) {
        log << strformat("Average %d profiles", profiles.size()) << std::endl;
    }
    if (distanceMatrix) {
        setCodeDist(out);
    }

    /* Compute constraints */
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < constraintSeqs.size(); i++) {
        for (size_t in = 0; in < profiles.size(); in++) {
            out.nOn[i] += profiles[in].nOn[i];
            out.nOff[i] += profiles[in].nOff[i];
        }
    }
}

AbsNeighbourJoining(void)::addToFreq(numeric_t fOut[], double weight, size_t codeIn, numeric_t fIn[]) {
    assert(fOut != NULL);
    if (fIn != nullptr) {
        operations.vector_add_mult(fOut, fIn, weight, options.nCodes);
    } else if (distanceMatrix) {
        assert(codeIn != NOCODE);
        operations.vector_add_mult(fOut, distanceMatrix.codeFreq[codeIn], weight, options.nCodes);
    } else {
        assert(codeIn != NOCODE);
        fOut[codeIn] += weight;
    }
}

/* Make the (unrotated) frequencies sum to 1
   Simply dividing by total_weight is not ideal because of roundoff error
   So compute total_freq instead
*/
AbsNeighbourJoining(void)::normalizeFreq(numeric_t freq[]) {
    double total_freq = 0;
    if (distanceMatrix) {
        /* The total frequency is dot_product(true_frequencies, 1)
           So we rotate the 1 vector by eigeninv (stored in eigentot)
        */
        total_freq = operations.vector_multiply_sum(freq, distanceMatrix.eigentot, options.nCodes);
    } else {
        for (int k = 0; k < options.nCodes; k++)
            total_freq += freq[k];
    }
    if (total_freq > options.fPostTotalTolerance) {
        numeric_t inverse_weight = 1.0 / total_freq;
        operations.vector_multiply_by(freq, inverse_weight, options.nCodes);
    } else {
        /* This can happen if we are in a very low-weight region, e.g. if a mostly-gap position gets weighted down
           repeatedly; just set them all to arbitrary but legal values */
        if (!distanceMatrix) {
            for (int k = 0; k < options.nCodes; k++)
                freq[k] = 1.0 / options.nCodes;
        } else {
            for (int k = 0; k < options.nCodes; k++)
                freq[k] = distanceMatrix.codeFreq[0][k];
        }
    }
}

AbsNeighbourJoining(void)::setCodeDist(Profile &profile) {
    if (profile.codeDist.empty()) {
        profile.codeDist.resize(nPos * options.nCodes);
    }
    size_t iFreq = 0;
    for (size_t i = 0; i < nPos; i++) {
        numeric_t *f = getFreq(profile, i, iFreq);

        for (int k = 0; k < options.nCodes; k++)
            profile.codeDist[i * options.nCodes + k] = profileDistPiece(profile.codes[i], k, f, NULL, NULL);
    }
    assert(iFreq == profile.vectors.size());
}

AbsNeighbourJoining(double)::
profileDistPiece(size_t code1, size_t code2, numeric_t f1[], numeric_t f2[], numeric_t codeDist2[]) {
    if (distanceMatrix) {
        if (code1 != NOCODE && code2 != NOCODE) { /* code1 vs code2 */
            return distanceMatrix.distances[code1][code2];
        } else if (codeDist2 != NULL && code1 != NOCODE) { /* code1 vs. codeDist2 */
            return codeDist2[code1];
        } else { /* f1 vs f2 */
            if (f1 == NULL) {
                if (code1 == NOCODE) { return (10.0); }
                f1 = &distanceMatrix.codeFreq[code1][0];
            }
            if (f2 == NULL) {
                if (code2 == NOCODE) { return (10.0); }
                f2 = &distanceMatrix.codeFreq[code2][0];
            }
            return operations.vector_multiply3_sum(f1, f2, distanceMatrix.eigenval, options.nCodes);
        }
    } else {
        /* no matrix */
        if (code1 != NOCODE) {
            if (code2 != NOCODE) {
                return (code1 == code2 ? 0.0 : 1.0); /* code1 vs code2 */
            } else {
                if (f2 == NULL) { return (10.0); }
                return 1.0 - f2[code1]; /* code1 vs. f2 */
            }
        } else {
            if (code2 != NOCODE) {
                if (f1 == NULL) { return (10.0); }
                return 1.0 - f1[code2]; /* f1 vs code2 */
            } else { /* f1 vs. f2 */
                if (f1 == NULL || f2 == NULL) { return (10.0); }
                double piece = 1.0;
                for (int k = 0; k < options.nCodes; k++) {
                    piece -= f1[k] * f2[k];
                }
                return piece;
            }
        }
    }
}

AbsNeighbourJoining(void)::setOutDistance(size_t iNode, size_t nActive) {
    if (nOutDistActive[iNode] == nActive) {
        return;
    }

    /* May be called by InitNJ before we have parents */
    assert(iNode >= 0 && (parent.empty() || parent[iNode] < 0));
    Besthit dist;
    profileDist(profiles[iNode], outprofile, dist);
    options.debug.outprofileOps++;

    /* out(A) = sum(X!=A) d(A,X)
       = sum(X!=A) (profiledist(A,X) - diam(A) - diam(X))
       = sum(X!=A) profiledist(A,X) - (N-1)*diam(A) - (totdiam - diam(A))

       in the absence of gaps:
       profiledist(A,out) = mean profiledist(A, all active nodes)
       sum(X!=A) profiledist(A,X) = N * profiledist(A,out) - profiledist(A,A)

       With gaps, we need to take the weights of the comparisons into account, where
       w(Ai) is the weight of position i in profile A:
       w(A,B) = sum_i w(Ai) * w(Bi)
       d(A,B) = sum_i w(Ai) * w(Bi) * d(Ai,Bi) / w(A,B)

       sum(X!=A) profiledist(A,X) ~= (N-1) * profiledist(A, Out w/o A)
       profiledist(A, Out w/o A) = sum_X!=A sum_i d(Ai,Xi) * w(Ai) * w(Bi) / ( sum_X!=A sum_i w(Ai) * w(Bi) )
       d(A, Out) = sum_A sum_i d(Ai,Xi) * w(Ai) * w(Bi) / ( sum_X sum_i w(Ai) * w(Bi) )

       and so we get
       profiledist(A,out w/o A) = (top of d(A,Out) - top of d(A,A)) / (weight of d(A,Out) - weight of d(A,A))
       top = dist * weight
       with another correction of nActive because the weight of the out-profile is the average
       weight not the total weight.
    */
    double top = (nActive - 1) * (dist.dist * dist.weight * nActive - selfweight[iNode] * selfdist[iNode]);
    double bottom = (dist.weight * nActive - selfweight[iNode]);
    double pdistOutWithoutA = top / bottom;
    outDistances[iNode] = bottom > 0.01 ?
                          pdistOutWithoutA - diameter[iNode] * (nActive - 1) -
                          (totdiam - diameter[iNode])
                                        : 3.0;
    nOutDistActive[iNode] = nActive;

    if (options.verbose > 3 && iNode < 5) {
        #pragma omp critical
        {
            log << strformat("NewOutDist for %d %f from dist %f selfd %f diam %f totdiam %f newActive %d",
                             iNode, outDistances[iNode], dist.dist, selfdist[iNode], diameter[iNode],
                             totdiam, nActive) << std::endl;
        };
    }
    if (options.verbose > 6 && (iNode % 10) == 0) {
        #pragma omp critical
        { ;
            /* Compute the actual out-distance and compare */
            double total = 0.0;
            double total_pd = 0.0;
            for (size_t j = 0; j < maxnode; j++) {
                if (j != iNode && (parent.empty() || parent[j] < 0)) {
                    Besthit bh;
                    profileDist(profiles[iNode], profiles[j], bh);
                    total_pd += bh.dist;
                    total += bh.dist - (diameter[iNode] + diameter[j]);
                }
            }
            log << strformat("OutDist for Node %d %f truth %f profiled %f truth %f pd_err %f",
                             iNode, outDistances[iNode], total, pdistOutWithoutA, total_pd,
                             fabs(pdistOutWithoutA - total_pd)) << std::endl;
        };
    }
}

AbsNeighbourJoining(void)::profileDist(Profile &profile1, Profile &profile2, Besthit &hit) {
    double top = 0;
    double denom = 0;
    size_t iFreq1 = 0;
    size_t iFreq2 = 0;
    for (size_t i = 0; i < nPos; i++) {
        numeric_t *f1 = getFreq(profile1, i, iFreq1);
        numeric_t *f2 = getFreq(profile2, i, iFreq2);
        if (profile1.weights[i] > 0 && profile2.weights[i] > 0) {
            double weight = profile1.weights[i] * profile2.weights[i];
            denom += weight;
            double piece = profileDistPiece(profile1.codes[i], profile2.codes[i], f1, f2,
                                            (!profile2.codeDist.empty() ? &profile2.codeDist[i * options.nCodes]
                                                                        : NULL));
            top += weight * piece;
        }
    }
    assert(iFreq1 == profile1.vectors.size());
    assert(iFreq2 == profile2.vectors.size());
    hit.weight = denom > 0 ? denom : 0.01; /* 0.01 is an arbitrarily low value of weight (normally >>1) */
    hit.dist = denom > 0 ? top / denom : 1;
    options.debug.profileOps++;
}

AbsNeighbourJoining(void)::seqDist(std::string &codes1, std::string &codes2, Besthit &hit) {
    double top = 0;        /* summed over positions */
    size_t nUse = 0;
    if (distanceMatrix) {
        int nDiff = 0;
        for (size_t i = 0; i < nPos; i++) {
            if (codes1[i] != NOCODE && codes2[i] != NOCODE) {
                nUse++;
                if (codes1[i] != codes2[i]) nDiff++;
            }
        }
        top = (double) nDiff;
    } else {
        for (size_t i = 0; i < nPos; i++) {
            if (codes1[i] != NOCODE && codes2[i] != NOCODE) {
                nUse++;
                top += distanceMatrix.distances[(size_t) codes1[i]][(size_t) codes2[i]];
            }
        }
    }
    hit.weight = (double) nUse;
    hit.dist = nUse > 0 ? top / (double) nUse : 1.0;
    options.debug.seqOps++;
}


AbsNeighbourJoining(inline Precision*)::getFreq(Profile &p, size_t &i, size_t &ivector) {
    return p.weights[i] > 0 && p.codes[i] == NOCODE ? &p.vectors[options.nCodes * (ivector++)] : NULL;
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