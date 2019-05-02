
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
            auto &constraintSeq = constraintSeqs[i];
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

AbsNeighbourJoining(void)::setProfile(size_t node, double weight1) {
    Children &c = child[node];
    assert(c.nChild == 2);
    assert(profiles.size() > c.child[0]);
    assert(profiles.size() > c.child[1]);
    averageProfile(profiles[node], profiles[c.child[0]], profiles[c.child[1]], weight1);
}

AbsNeighbourJoining(void)::averageProfile(Profile &out, Profile &profile1, Profile &profile2, double bionjWeight) {
    if (bionjWeight < 0) {
        bionjWeight = 0.5;
    }
    if (!out.codeDist.empty()) {
        out.codeDist.clear();
    }

    for (size_t i = 0; i < nPos; i++) {
        out.weights[i] = bionjWeight * profile1.weights[i] + (1 - bionjWeight) * profile2.weights[i];
        out.codes[i] = NOCODE;
        if (out.weights[i] > 0) {
            if (profile1.weights[i] > 0 && profile1.codes[i] != NOCODE
                && (profile2.weights[i] <= 0 || profile1.codes[i] == profile2.codes[i])) {
                out.codes[i] = profile1.codes[i];
            } else if (profile1.weights[i] <= 0
                       && profile2.weights[i] > 0
                       && profile2.codes[i] != NOCODE) {
                out.codes[i] = profile2.codes[i];
            }
            if (out.codes[i] == NOCODE) {
                out.vectors.resize(out.vectors.size() + 1);
            }
        }
    }

    /* Allocate and set the vectors */
    assert(out.vectors.size() == options.nCodes * out.vectors.size());
    for (size_t i = 0; i < options.nCodes * out.vectors.size(); i++) {
        out.vectors[i] = 0;
    }
    options.debug.nProfileFreqAlloc += out.vectors.size();
    options.debug.nProfileFreqAvoid += nPos - out.vectors.size();
    size_t iFreqOut = 0;
    size_t iFreq1 = 0;
    size_t iFreq2 = 0;
    for (size_t i = 0; i < nPos; i++) {
        numeric_t *f = getFreq(out, i, iFreqOut);
        numeric_t *f1 = getFreq(profile1, i, iFreq1);
        numeric_t *f2 = getFreq(profile2, i, iFreq2);
        if (f != nullptr) {
            if (profile1.weights[i] > 0) {
                addToFreq(f, profile1.weights[i] * bionjWeight, profile1.codes[i], f1);
            }
            if (profile2.weights[i] > 0) {
                addToFreq(f, profile2.weights[i] * (1.0 - bionjWeight), profile2.codes[i], f2);
            }
            normalizeFreq(f);
        } /* end if computing f */
        if (options.verbose > 10 && i < 5) {
            log << strformat("Average profiles: pos %d in-w1 %f in-w2 %f bionjWeight %f to weight %f code %d",
                             i, profile1.weights[i], profile2.weights[i], bionjWeight, out.weights[i], out.codes[i]);
            if (f != NULL) {
                for (int k = 0; k < options.nCodes; k++) {
                    log << strformat("\t%c:%f", options.codesString[k], f ? f[k] : -1.0);
                }
                log << std::endl;
            }
        }
    } /* end loop over positions */
    assert(iFreq1 == profile1.vectors.size());
    assert(iFreq2 == profile2.vectors.size());
    assert(iFreqOut == out.vectors.size());

    /* compute total constraints */
    for (size_t i = 0; i < constraintSeqs.size(); i++) {
        out.nOn[i] = profile1.nOn[i] + profile2.nOn[i];
        out.nOff[i] = profile1.nOff[i] + profile2.nOff[i];
    }
    options.debug.profileAvgOps++;
}

AbsNeighbourJoining(void)::readTree(Uniquify &unique, HashTable &hashnames, std::istream &fpInTree) {
    assert(seqs.size() == unique.uniqueSeq.size());
    /* First, do a preliminary parse of the tree to with non-unique leaves ignored
       We need to store this separately from NJ because it may have too many internal nodes
       (matching sequences show up once in the NJ but could be in multiple places in the tree)
       Will use iUnique as the index of nodes, as in the NJ structure
    */
    assert(maxnodes == unique.alnToUniq.size() * 2);
    assert(maxnode == unique.alnToUniq.size());
    std::vector<int64_t> parents(maxnodes, -1);
    std::vector<Children> children(maxnodes);

    /* The stack is the current path to the root, with the root at the first (top) position */
    size_t stack_size = 1;
    std::vector<size_t> stack(maxnodes);
    stack[0] = root;
    size_t nDown = 0;
    size_t nUp = 0;

    std::string token;
    token.reserve(5000);

    if (readTreeToken(fpInTree, token) || token[0] != '(') {
        readTreeError("No '(' at start", token);
    }
    /* nDown is still 0 because we have created the root */

    while (readTreeToken(fpInTree, token)) {
        if (nDown > 0) {        /* In a stream of parentheses */
            if (token[0] == '(') {
                nDown++;
            } else if (token[0] == ',' || token[0] == ';' || token[0] == ':' || token[0] == ')') {
                readTreeError("while reading parentheses", token);
            } else {
                /* Add intermediate nodes if nDown was > 1 (for nDown=1, the only new node is the leaf) */
                while (nDown-- > 0) {
                    size_t newnode = maxnode++;
                    assert(newnode < maxnodes);
                    readTreeAddChild(stack[stack_size - 1], newnode, parents, children);
                    if (options.verbose > 5) {
                        log << strformat("Added internal child %d of %d, stack size increase to %d",
                                         newnode, stack[stack_size - 1], stack_size + 1) << std::endl;
                    }
                    stack[stack_size++] = newnode;
                    assert(stack_size < maxnodes);
                }
                readTreeMaybeAddLeaf(stack[stack_size - 1], token, hashnames, unique, parents, children);
            }
        } else if (nUp > 0) {
            if (token[0] == ';') {    /* end the tree? */
                if (nUp != stack_size) {
                    readTreeError("unbalanced parentheses", token);
                } else {
                    break;
                }
            } else if (token[0] == ')') {
                nUp++;
            } else if (token[0] == '(') {
                readTreeError("unexpected '(' after ')'", token);
            } else if (token[0] == ':') {
                /* Read the branch length and ignore it */
                if (!readTreeToken(fpInTree, token) || (token[0] != '-' && !std::isdigit(token[0]))) {
                    readTreeError("not recognized as a branch length", token);
                }
            } else if (token[0] == ',') {
                /* Go back up the stack the correct #times */
                while (nUp-- > 0) {
                    stack_size--;
                    if (options.verbose > 5) {
                        log << strformat("Up to nUp=%d stack size %d at %d",
                                         nUp, stack_size, stack[stack_size - 1]) << std::endl;
                    }
                    if (stack_size <= 0) {
                        readTreeError("too many ')'", token);
                    }
                }
                nUp = 0;
            } else if (token[0] == '-' || std::isdigit(token[0])) {
                /* ignore bootstrap value */
            } else {
                log << "Warning while parsing tree: non-numeric label " << token << " for internal node" << std::endl;
            }
        } else if (token[0] == '(') {
            nDown = 1;
        } else if (token[0] == ')') {
            nUp = 1;
        } else if (token[0] == ':') {
            if (!readTreeToken(fpInTree, token) || (token[0] != '-' && !std::isdigit(token[0]))) {
                readTreeError("not recognized as a branch length", token);
            }
        } else if (token[0] == ',') {
            /* do nothing */
        } else if (token[0] == ';') {
            readTreeError("unexpected token", token);
        } else {
            readTreeMaybeAddLeaf(stack[stack_size - 1], token, hashnames, unique, parents, children);
        }
    }

    /* Verify that all sequences were seen */
    for (size_t i = 0; i < unique.uniqueSeq.size(); i++) {
        if (parents[i] < 0) {
            throw std::invalid_argument(
                    strformat("Alignment sequence %d (unique %d) absent from input tree\n"
                              "The starting tree (the argument to -intree) must include all sequences in the alignment!",
                              unique.uniqueFirst[i], i));
        }
    }

    /* Simplify the tree -- remove all internal nodes with < 2 children
       Keep trying until no nodes get removed
    */
    size_t nRemoved;
    do {
        nRemoved = 0;
        /* Here stack is the list of nodes we haven't visited yet while doing
           a tree traversal */
        stack_size = 1;
        stack[0] = root;
        while (stack_size > 0) {
            int node = stack[--stack_size];
            if (node >= (int64_t) unique.uniqueSeq.size()) { /* internal node */
                if (children[node].nChild <= 1) {
                    if (node != root) {
                        readTreeRemove(parents, children, node);
                        nRemoved++;
                    } else if (node == root && children[node].nChild == 1) {
                        size_t newroot = children[node].child[0];
                        parents[newroot] = -1;
                        children[root].nChild = 0;
                        nRemoved++;
                        if (options.verbose > 5) {
                            log << strformat("Changed root from %d to %d", root, newroot) << std::endl;
                        }
                        root = newroot;
                        stack[stack_size++] = newroot;
                    }
                } else {
                    int j;
                    for (j = 0; j < children[node].nChild; j++) {
                        assert(stack_size < maxnodes);
                        stack[stack_size++] = children[node].child[j];
                        if (options.verbose > 5) {
                            log << strformat("Added %d to stack", stack[stack_size - 1]) << std::endl;
                        }
                    }
                }
            }
        }
    } while (nRemoved > 0);

    /* Simplify the root node to 3 children if it has 2 */
    if (children[root].nChild == 2) {
        for (size_t i = 0; i < 2; i++) {
            size_t child = children[root].child[i];
            assert(child >= 0 && child < maxnodes);
            if (children[child].nChild == 2) {
                readTreeRemove(parents, children, child); /* replace root -> child -> A,B with root->A,B */
                break;
            }
        }
    }

    for (size_t i = 0; i < maxnodes; i++)
        if (options.verbose > 5) {
            log << strformat("Simplfied node %d has parent %d nchild %d", i, parents[i], children[i].nChild)
                << std::endl;
        }

    /* Map the remaining internal nodes to NJ nodes */
    std::vector<size_t> map(maxnodes);
    for (size_t i = 0; i < unique.uniqueSeq.size(); i++) {
        map[i] = i;
    }
    for (size_t i = unique.uniqueSeq.size(); i < maxnodes; i++) {
        map[i] = -1;
    }
    stack_size = 1;
    stack[0] = root;
    while (stack_size > 0) {
        int64_t node = stack[--stack_size];
        if (node >= (int64_t) unique.uniqueSeq.size()) { /* internal node */
            assert(node == root || children[node].nChild > 1);
            map[node] = maxnode++;
            for (int64_t i = 0; i < children[node].nChild; i++) {
                assert(stack_size < maxnodes);
                stack[stack_size++] = children[node].child[i];
            }
        }
    }
    for (size_t i = 0; i < maxnodes; i++)
        if (options.verbose > 5) {
            log << strformat("Map %d to %d (parent %d nchild %d)", i, map[i], parents[i], children[i].nChild)
                << std::endl;
        }

    /* Set parents, children, root */
    root = map[root];
    int64_t node;
    for (node = 0; node < (int64_t) maxnodes; node++) {
        int njnode = map[node];
        if (njnode >= 0) {
            child[njnode].nChild = children[node].nChild;
            for (int64_t i = 0; i < children[node].nChild; i++) {
                assert(children[node].child[i] >= 0 && children[node].child[i] < maxnodes);
                child[njnode].child[i] = map[children[node].child[i]];
            }
            if (parents[node] >= 0) {
                this->parent[njnode] = map[parents[node]];
            }
        }
    }

    /* Make sure that parents/child relationships match */
    for (size_t i = 0; i < maxnode; i++) {
        Children &c = child[i];
        for (int64_t j = 0; j < c.nChild; j++) {
            assert(c.child[j] >= 0 && c.child[j] < maxnode && this->parent[c.child[j]] == (int64_t) i);
        }
    }
    assert(this->parent[root] < 0);

    /* Compute profiles as balanced -- the NNI stage will recompute these
       profiles anyway
    */
    std::vector<bool> traversal(maxnode, false);
    node = root;
    while ((node = TraversePostorder(node, traversal, nullptr)) >= 0) {
        if (node >= (int64_t) seqs.size() && node != root) {
            setProfile(node, -1.0);
        }
    }
}

AbsNeighbourJoining(void)::
printNJ(std::ostream &out, std::vector<std::string> &names, Uniquify &unique, bool bShowSupport) {
    /* And print the tree: depth first search
     * The stack contains
     * list of remaining children with their depth
     * parent node, with a flag of -1 so I know to print right-paren
     */
    constexpr auto FP_FORMAT = sizeof(numeric_t)==sizeof(float)? "%.5f" : "%.9f";

    if (seqs.size() == 1 && unique.alnNext[unique.uniqueFirst[0]] >= 0) {
        /* Special case -- otherwise we end up with double parens */
        size_t first = unique.uniqueFirst[0];
        assert(first >= 0 && first < seqs.size());
        out << "(";
        quotes(out, names[first], options.bQuote) << ":0.0";
        int64_t iName = unique.alnNext[first];
        while (iName >= 0) {
            assert(iName < (int64_t)seqs.size());
            out << ",";
            quotes(out, names[iName], options.bQuote) << ":0.0";
            iName = unique.alnNext[iName];
        }
        out << ");" << std::endl;
        return;
    }

    size_t stackSize = 1;
    std::vector<std::pair<int64_t, int64_t>> stack(maxnodes);
    stack[0].first = root; //node
    stack[0].second = 0; //end

    while (stackSize > 0) {
        std::pair<int64_t, int64_t> &last = stack[stackSize - 1];
        stackSize--;
        /* Save last, as we are about to overwrite it */
        int64_t node = last.first;
        int64_t end = last.second;

        if (node < (int64_t)seqs.size()) {
            if ((int64_t)child[parent[node]].child[0] != node) {
                out << ",";
            }
            int64_t first = unique.uniqueFirst[node];
            assert(first >= 0 && first < (int64_t)seqs.size());
            /* Print the name, or the subtree of duplicate names */
            if (unique.alnNext[first] == -1) {
                quotes(out, names[first], options.bQuote);
            } else {
                out << "(";
                quotes(out, names[first], options.bQuote) << ":0.0";
                int64_t iName = unique.alnNext[first];
                while (iName >= 0) {
                    assert(iName < (int64_t)seqs.size());
                    out << ",";
                    quotes(out, names[iName], options.bQuote) << ":0.0";
                    iName = unique.alnNext[iName];
                }
                out << ")";
            }

            /* Print the branch length */
            out << ":" << strformat(FP_FORMAT, branchlength[node]);
        } else if (end) {
            if (node == root) {
                out << ")";
            }else if (bShowSupport) {
                out << strformat(")%.3f:",support[node]) << strformat(FP_FORMAT, branchlength[node]);
            }else{
                out << "):" << strformat(FP_FORMAT, branchlength[node]);
            }
        } else {
            if (node != root && (int64_t)child[parent[node]].child[0] != node){
                out << ",";
            }
            out << "(";
            stackSize++;
            stack[stackSize - 1].first = node;
            stack[stackSize - 1].second = 1;
            Children &c = child[node];
            /* put children on in reverse order because we use the last one first */
            int i;
            for (i = c.nChild - 1; i >= 0; i--) {
                stackSize++;
                stack[stackSize - 1].first = c.child[i];
                stack[stackSize - 1].second = 0;
            }
        }
    }
    out << ";" << std::endl;
}

AbsNeighbourJoining(void)::fastNJ() {
    #if 0 //TODO translate
    int64_t iNode;

    assert(seqs.size() >= 1);
    if (seqs.size() < 3) {
        root = maxnode++;
        child[root].nChild = seqs.size();
        for (iNode = 0; iNode < seqs.size(); iNode++) {
            parent[iNode] = root;
            child[root].child[iNode] = iNode;
        }
        if (seqs.size() == 1) {
            branchlength[0] = 0;
        } else {
            assert (seqs.size() == 2);
            Besthit hit;
            seqDist(profiles[0]->codes,profiles[1]->codes,hit);
            branchlength[0] = hit.dist/2.0;
            branchlength[1] = hit.dist/2.0;
        }
        return;
    }

    /* else 3 or more sequences */

    /* The visible set stores the best hit of each node (unless using top hits, in which case
       it is handled by the top hits routines) */
    besthit_t *visible = NULL;	/* Not used if doing top hits */
    besthit_t *besthitNew = NULL;	/* All hits of new node -- not used if doing top-hits */

    /* The top-hits lists, with the key parameter m = length of each top-hit list */
    top_hits_t *tophits = NULL;
    int m = 0;			/* maximum length of a top-hits list */
    if (tophitsMult > 0) {
        m = (int)(0.5 + tophitsMult*sqrt(seqs.size()));
        if(m<4 || 2*m >= seqs.size()) {
            m=0;
            if(options.verbose>1) fprintf(stderr,"Too few leaves, turning off top-hits\n");
        } else {
            if(options.verbose>2) fprintf(stderr,"Top-hit-list size = %d of %d\n", m, seqs.size());
        }
    }
    assert(!(options.slow && m>0));

    /* Initialize top-hits or visible set */
    if (m>0) {
        tophits = InitTopHits(NJ, m);
        SetAllLeafTopHits(/*IN/UPDATE*/NJ, /*OUT*/tophits);
        ResetTopVisible(/*IN/UPDATE*/NJ, /*nActive*/seqs.size(), /*IN/OUT*/tophits);
    } else if (!slow) {
        visible = (besthit_t*)mymalloc(sizeof(besthit_t)*NJ->maxnodes);
        besthitNew = (besthit_t*)mymalloc(sizeof(besthit_t)*NJ->maxnodes);
        for (iNode = 0; iNode < seqs.size(); iNode++)
            SetBestHit(iNode, NJ, /*nActive*/seqs.size(), /*OUT*/&visible[iNode], /*OUT IGNORED*/NULL);
    }

    /* Iterate over joins */
    int nActiveOutProfileReset = seqs.size();
    int nActive;
    for (nActive = seqs.size(); nActive > 3; nActive--) {
        int nJoinsDone = seqs.size() - nActive;
        if (nJoinsDone > 0 && (nJoinsDone % 100) == 0)
            ProgressReport("Joined %6d of %6d", nJoinsDone, seqs.size()-3, 0, 0);

        besthit_t join; 		/* the join to do */
        if (slow) {
            ExhaustiveNJSearch(NJ,nActive,/*OUT*/&join);
        } else if (m>0) {
            TopHitNJSearch(/*IN/UPDATE*/NJ, nActive, /*IN/OUT*/tophits, /*OUT*/&join);
        } else {
            FastNJSearch(NJ, nActive, /*IN/OUT*/visible, /*OUT*/&join);
        }

        if (options.verbose>2) {
            double penalty = constraintWeight
                             * (double)JoinConstraintPenalty(NJ, join.i, join.j);
            if (penalty > 0.001) {
                fprintf(stderr, "Constraint violation during neighbor-joining %d %d into %d penalty %.3f\n",
                        join.i, join.j, NJ->maxnode, penalty);
                int iC;
                for (iC = 0; iC < NJ->nConstraints; iC++) {
                    int local = JoinConstraintPenaltyPiece(NJ, join.i, join.j, iC);
                    if (local > 0)
                        fprintf(stderr, "Constraint %d piece %d %d/%d %d/%d %d/%d\n", iC, local,
                                NJ->profiles[join.i]->nOn[iC],
                                NJ->profiles[join.i]->nOff[iC],
                                NJ->profiles[join.j]->nOn[iC],
                                NJ->profiles[join.j]->nOff[iC],
                                NJ->outprofile->nOn[iC] - NJ->profiles[join.i]->nOn[iC] - NJ->profiles[join.j]->nOn[iC],
                                NJ->outprofile->nOff[iC] - NJ->profiles[join.i]->nOff[iC] - NJ->profiles[join.j]->nOff[iC]);
                }
            }
        }

        /* because of the stale out-distance heuristic, make sure that these are up-to-date */
        SetOutDistance(NJ, join.i, nActive);
        SetOutDistance(NJ, join.j, nActive);
        /* Make sure weight is set and criterion is up to date */
        SetDistCriterion(NJ, nActive, /*IN/OUT*/&join);
        assert(NJ->nOutDistActive[join.i] == nActive);
        assert(NJ->nOutDistActive[join.j] == nActive);

        int newnode = NJ->maxnode++;
        parent[join.i] = newnode;
        parent[join.j] = newnode;
        child[newnode].nChild = 2;
        child[newnode].child[0] = join.i < join.j ? join.i : join.j;
        child[newnode].child[1] = join.i > join.j ? join.i : join.j;

        double rawIJ = join.dist + NJ->diameter[join.i] + NJ->diameter[join.j];
        double distIJ = join.dist;

        double deltaDist = (NJ->outDistances[join.i]-NJ->outDistances[join.j])/(double)(nActive-2);
        NJ->branchlength[join.i] = (distIJ + deltaDist)/2;
        NJ->branchlength[join.j] = (distIJ - deltaDist)/2;

        double bionjWeight = 0.5;	/* IJ = bionjWeight*I + (1-bionjWeight)*J */
        double varIJ = rawIJ - NJ->varDiameter[join.i] - NJ->varDiameter[join.j];

        if (bionj && join.weight > 0.01 && varIJ > 0.001) {
            /* Set bionjWeight according to the BIONJ formula, where
           the variance matrix is approximated by
      
           Vij = ProfileVar(i,j) - varDiameter(i) - varDiameter(j)
           ProfileVar(i,j) = distance(i,j) = top(i,j)/weight(i,j)
      
           (The node's distance diameter does not affect the variances.)
      
           The BIONJ formula is equation 9 from Gascuel 1997:
      
           bionjWeight = 1/2 + sum(k!=i,j) (Vjk - Vik) / ((nActive-2)*Vij)
           sum(k!=i,j) (Vjk - Vik) = sum(k!=i,j) Vik - varDiameter(j) + varDiameter(i)
           = sum(k!=i,j) ProfileVar(j,k) - sum(k!=i,j) ProfileVar(i,k) + (nActive-2)*(varDiameter(i)-varDiameter(j))
      
           sum(k!=i,j) ProfileVar(i,k)
           ~= (sum(k!=i,j) distance(i,k) * weight(i,k))/(mean(k!=i,j) weight(i,k))
           ~= (N-2) * top(i, Out-i-j) / weight(i, Out-i-j)
      
           weight(i, Out-i-j) = N*weight(i,Out) - weight(i,i) - weight(i,j)
           top(i, Out-i-j) = N*top(i,Out) - top(i,i) - top(i,j)
            */
            besthit_t outI;
            besthit_t outJ;
            ProfileDist(NJ->profiles[join.i],NJ->outprofile,NJ->nPos,NJ->distance_matrix,/*OUT*/&outI);
            ProfileDist(NJ->profiles[join.j],NJ->outprofile,NJ->nPos,NJ->distance_matrix,/*OUT*/&outJ);
            outprofileOps += 2;

            double varIWeight = (nActive * outI.weight - NJ->selfweight[join.i] - join.weight);
            double varJWeight = (nActive * outJ.weight - NJ->selfweight[join.j] - join.weight);

            double varITop = outI.dist * outI.weight * nActive
                             - NJ->selfdist[join.i] * NJ->selfweight[join.i] - rawIJ * join.weight;
            double varJTop = outJ.dist * outJ.weight * nActive
                             - NJ->selfdist[join.j] * NJ->selfweight[join.j] - rawIJ * join.weight;

            double deltaProfileVarOut = (nActive-2) * (varJTop/varJWeight - varITop/varIWeight);
            double deltaVarDiam = (nActive-2)*(NJ->varDiameter[join.i] - NJ->varDiameter[join.j]);
            if (varJWeight > 0.01 && varIWeight > 0.01)
                bionjWeight = 0.5 + (deltaProfileVarOut+deltaVarDiam)/(2*(nActive-2)*varIJ);
            if(bionjWeight<0) bionjWeight=0;
            if(bionjWeight>1) bionjWeight=1;
            if (verbose>2) fprintf(stderr,"dVarO %f dVarDiam %f varIJ %f from dist %f weight %f (pos %d) bionjWeight %f %f\n",
                                   deltaProfileVarOut, deltaVarDiam,
                                   varIJ, join.dist, join.weight, NJ->nPos,
                                   bionjWeight, 1-bionjWeight);
            if (verbose>3 && (newnode%5) == 0) {
                /* Compare weight estimated from outprofiles from weight made by summing over other nodes */
                double deltaProfileVarTot = 0;
                for (iNode = 0; iNode < newnode; iNode++) {
                    if (parent[iNode] < 0) { /* excludes join.i, join.j */
                        besthit_t di, dj;
                        ProfileDist(NJ->profiles[join.i],NJ->profiles[iNode],NJ->nPos,NJ->distance_matrix,/*OUT*/&di);
                        ProfileDist(NJ->profiles[join.j],NJ->profiles[iNode],NJ->nPos,NJ->distance_matrix,/*OUT*/&dj);
                        deltaProfileVarTot += dj.dist - di.dist;
                    }
                }
                double lambdaTot = 0.5 + (deltaProfileVarTot+deltaVarDiam)/(2*(nActive-2)*varIJ);
                if (lambdaTot < 0) lambdaTot = 0;
                if (lambdaTot > 1) lambdaTot = 1;
                if (fabs(bionjWeight-lambdaTot) > 0.01 || verbose > 4)
                    fprintf(stderr, "deltaProfileVar actual %.6f estimated %.6f lambda actual %.3f estimated %.3f\n",
                            deltaProfileVarTot,deltaProfileVarOut,lambdaTot,bionjWeight);
            }
        }
        if (verbose > 2) fprintf(stderr, "Join\t%d\t%d\t%.6f\tlambda\t%.6f\tselfw\t%.3f\t%.3f\tnew\t%d\n",
                                 join.i < join.j ? join.i : join.j,
                                 join.i < join.j ? join.j : join.i,
                                 join.criterion, bionjWeight,
                                 NJ->selfweight[join.i < join.j ? join.i : join.j],
                                 NJ->selfweight[join.i < join.j ? join.j : join.i],
                                 newnode);

        NJ->diameter[newnode] = bionjWeight * (NJ->branchlength[join.i] + NJ->diameter[join.i])
                                + (1-bionjWeight) * (NJ->branchlength[join.j] + NJ->diameter[join.j]);
        NJ->varDiameter[newnode] = bionjWeight * NJ->varDiameter[join.i]
                                   + (1-bionjWeight) * NJ->varDiameter[join.j]
                                   + bionjWeight * (1-bionjWeight) * varIJ;

        NJ->profiles[newnode] = AverageProfile(NJ->profiles[join.i],NJ->profiles[join.j],
                                               NJ->nPos, NJ->nConstraints,
                                               NJ->distance_matrix,
                                               bionj ? bionjWeight : /*noweight*/-1.0);

        /* Update out-distances and total diameters */
        int changedActiveOutProfile = nActiveOutProfileReset - (nActive-1);
        if (changedActiveOutProfile >= nResetOutProfile
            && changedActiveOutProfile >= fResetOutProfile * nActiveOutProfileReset) {
            /* Recompute the outprofile from scratch to avoid roundoff error */
            profile_t **activeProfiles = (profile_t**)mymalloc(sizeof(profile_t*)*(nActive-1));
            int nSaved = 0;
            NJ->totdiam = 0;
            for (iNode=0;iNode<NJ->maxnode;iNode++) {
                if (parent[iNode]<0) {
                    assert(nSaved < nActive-1);
                    activeProfiles[nSaved++] = NJ->profiles[iNode];
                    NJ->totdiam += NJ->diameter[iNode];
                }
            }
            assert(nSaved==nActive-1);
            FreeProfile(NJ->outprofile, NJ->nPos, NJ->nConstraints);
            if(verbose>2) fprintf(stderr,"Recomputing outprofile %d %d\n",nActiveOutProfileReset,nActive-1);
            NJ->outprofile = OutProfile(activeProfiles, nSaved,
                                        NJ->nPos, NJ->nConstraints,
                                        NJ->distance_matrix);
            activeProfiles = myfree(activeProfiles, sizeof(profile_t*)*(nActive-1));
            nActiveOutProfileReset = nActive-1;
        } else {
            UpdateOutProfile(/*OUT*/NJ->outprofile,
                                    NJ->profiles[join.i], NJ->profiles[join.j], NJ->profiles[newnode],
                                    nActive,
                                    NJ->nPos, NJ->nConstraints,
                                    NJ->distance_matrix);
            NJ->totdiam += NJ->diameter[newnode] - NJ->diameter[join.i] - NJ->diameter[join.j];
        }

        /* Store self-dist for use in other computations */
        besthit_t selfdist;
        ProfileDist(NJ->profiles[newnode],NJ->profiles[newnode],NJ->nPos,NJ->distance_matrix,/*OUT*/&selfdist);
        NJ->selfdist[newnode] = selfdist.dist;
        NJ->selfweight[newnode] = selfdist.weight;

        /* Find the best hit of the joined node IJ */
        if (m>0) {
            TopHitJoin(newnode, /*IN/UPDATE*/NJ, nActive-1, /*IN/OUT*/tophits);
        } else {
            /* Not using top-hits, so we update all out-distances */
            for (iNode = 0; iNode < NJ->maxnode; iNode++) {
                if (parent[iNode] < 0) {
                    /* True nActive is now nActive-1 */
                    SetOutDistance(/*IN/UPDATE*/NJ, iNode, nActive-1);
                }
            }

            if(visible != NULL) {
                SetBestHit(newnode, NJ, nActive-1, /*OUT*/&visible[newnode], /*OUT OPTIONAL*/besthitNew);
                if (verbose>2)
                    fprintf(stderr,"Visible %d %d %f %f\n",
                            visible[newnode].i, visible[newnode].j,
                            visible[newnode].dist, visible[newnode].criterion);
                if (besthitNew != NULL) {
                    /* Use distances to new node to update visible set entries that are non-optimal */
                    for (iNode = 0; iNode < NJ->maxnode; iNode++) {
                        if (parent[iNode] >= 0 || iNode == newnode)
                            continue;
                        int iOldVisible = visible[iNode].j;
                        assert(iOldVisible>=0);
                        assert(visible[iNode].i == iNode);

                        /* Update the criterion; use nActive-1 because haven't decremented nActive yet */
                        if (parent[iOldVisible] < 0)
                            SetCriterion(/*IN/OUT*/NJ, nActive-1, &visible[iNode]);

                        if (parent[iOldVisible] >= 0
                            || besthitNew[iNode].criterion < visible[iNode].criterion) {
                            if(verbose>3) fprintf(stderr,"Visible %d reset from %d to %d (%f vs. %f)\n",
                                                  iNode, iOldVisible,
                                                  newnode, visible[iNode].criterion, besthitNew[iNode].criterion);
                            if(parent[iOldVisible] < 0) nVisibleUpdate++;
                            visible[iNode].j = newnode;
                            visible[iNode].dist = besthitNew[iNode].dist;
                            visible[iNode].criterion = besthitNew[iNode].criterion;
                        }
                    } /* end loop over all nodes */
                } /* end if recording all hits of new node */
            } /* end if keeping a visible set */
        } /* end else (m==0) */
    } /* end loop over nActive */

#ifdef TRACK_MEMORY
    if (verbose>1) {
    struct mallinfo mi = mallinfo();
    fprintf(stderr, "Memory @ end of FastNJ(): %.2f MB (%.1f byte/pos) useful %.2f expected %.2f\n",
	    (mi.arena+mi.hblkhd)/1.0e6, (mi.arena+mi.hblkhd)/(double)(seqs.size()*(double)NJ->nPos),
	    mi.uordblks/1.0e6, mymallocUsed/1e6);
  }
#endif

    /* We no longer need the tophits, visible set, etc. */
    if (visible != NULL) visible = myfree(visible,sizeof(besthit_t)*NJ->maxnodes);
    if (besthitNew != NULL) besthitNew = myfree(besthitNew,sizeof(besthit_t)*NJ->maxnodes);
    tophits = FreeTopHits(tophits);

    /* Add a root for the 3 remaining nodes */
    int top[3];
    int nTop = 0;
    for (iNode = 0; iNode < NJ->maxnode; iNode++) {
        if (parent[iNode] < 0) {
            assert(nTop <= 2);
            top[nTop++] = iNode;
        }
    }
    assert(nTop==3);

    root = NJ->maxnode++;
    child[root].nChild = 3;
    for (nTop = 0; nTop < 3; nTop++) {
        parent[top[nTop]] = root;
        child[root].child[nTop] = top[nTop];
    }

    besthit_t dist01, dist02, dist12;
    ProfileDist(NJ->profiles[top[0]], NJ->profiles[top[1]], NJ->nPos, NJ->distance_matrix, /*OUT*/&dist01);
    ProfileDist(NJ->profiles[top[0]], NJ->profiles[top[2]], NJ->nPos, NJ->distance_matrix, /*OUT*/&dist02);
    ProfileDist(NJ->profiles[top[1]], NJ->profiles[top[2]], NJ->nPos, NJ->distance_matrix, /*OUT*/&dist12);

    double d01 = dist01.dist - NJ->diameter[top[0]] - NJ->diameter[top[1]];
    double d02 = dist02.dist - NJ->diameter[top[0]] - NJ->diameter[top[2]];
    double d12 = dist12.dist - NJ->diameter[top[1]] - NJ->diameter[top[2]];
    NJ->branchlength[top[0]] = (d01 + d02 - d12)/2;
    NJ->branchlength[top[1]] = (d01 + d12 - d02)/2;
    NJ->branchlength[top[2]] = (d02 + d12 - d01)/2;

    /* Check how accurate the outprofile is */
    if (verbose>2) {
        profile_t *p[3] = {NJ->profiles[top[0]], NJ->profiles[top[1]], NJ->profiles[top[2]]};
        profile_t *out = OutProfile(p, 3, NJ->nPos, NJ->nConstraints, NJ->distance_matrix);
        int i;
        double freqerror = 0;
        double weighterror = 0;
        for (i=0;i<NJ->nPos;i++) {
            weighterror += fabs(out->weights[i] - NJ->outprofile->weights[i]);
            int k;
            for(k=0;k<nCodes;k++)
                freqerror += fabs(out->vectors[nCodes*i+k] - NJ->outprofile->vectors[nCodes*i+k]);
        }
        fprintf(stderr,"Roundoff error in outprofile@end: WeightError %f FreqError %f\n", weighterror, freqerror);
        FreeProfile(out, NJ->nPos, NJ->nConstraints);
    }
    return;
    #endif
}

AbsNeighbourJoining(void)::
readTreeAddChild(size_t parent, size_t child, std::vector<int64_t> &parents, std::vector<Children> &children) {
    assert(parent >= 0);
    assert(child >= 0);
    assert(parents[child] < 0);
    assert(children[parent].nChild < 3);
    parents[child] = parent;
    children[parent].child[children[parent].nChild++] = child;
}

AbsNeighbourJoining(void)::readTreeMaybeAddLeaf(size_t parent, std::string &name, HashTable &hashnames,
                                                Uniquify &unique,
                                                std::vector<int64_t> &parents, std::vector<Children> &children) {
    auto hi = hashnames.find(name);
    if (hi == nullptr)
        readTreeError("not recognized as a sequence name", name);

    size_t iSeqNonunique = *hi;
    assert(iSeqNonunique >= 0 && iSeqNonunique < unique.alnToUniq.size());
    size_t iSeqUnique = unique.alnToUniq[iSeqNonunique];
    assert(iSeqUnique >= 0 && iSeqUnique < unique.uniqueSeq.size());
    /* Either record this leaves' parent (if it is -1) or ignore this leaf (if already seen) */
    if (parents[iSeqUnique] < 0) {
        readTreeAddChild(parent, iSeqUnique, /*IN/OUT*/parents, /*IN/OUT*/children);
        if (options.verbose > 5) {
            log << strformat("Found leaf uniq%d name %s child of %d", iSeqUnique, name.c_str(), parent) << std::endl;
        }
    } else {
        if (options.verbose > 5) {
            log << strformat("Skipped redundant leaf uniq%d name %s", iSeqUnique, name.c_str()) << std::endl;
        }
    }
}

AbsNeighbourJoining(void)::readTreeRemove(std::vector<int64_t> &parents, std::vector<Children> &children, size_t node) {
    if (options.verbose > 5) {
        log << strformat("Removing node %d parent %d", node, parents[node]) << std::endl;
    }
    assert(parents[node] >= 0);
    int64_t parent = parents[node];
    parents[node] = -1;
    Children &pc = children[parent];
    int oldn = 0;
    for (; oldn < pc.nChild; oldn++) {
        if (pc.child[oldn] == node)
            break;
    }
    assert(oldn < pc.nChild);

    /* move successor nodes back in child list and shorten list */
    for (int i = oldn; i < pc.nChild - 1; i++) {
        pc.child[i] = pc.child[i + 1];
    }
    pc.nChild--;

    /* add its children to parent's child list */
    Children &nc = children[node];
    if (nc.nChild > 0) {
        assert(nc.nChild <= 2);
        assert(pc.nChild < 3);
        assert(pc.nChild + nc.nChild <= 3);
        for (int j = 0; j < nc.nChild; j++) {
            if (options.verbose > 5) {
                log << strformat("Repointing parent %d to child %d", parent, nc.child[j]) << std::endl;
            }
            pc.child[pc.nChild++] = nc.child[j];
            parents[nc.child[j]] = parent;
        }
        nc.nChild = 0;
    }
}

AbsNeighbourJoining(bool)::readTreeToken(std::istream &fpInTree, std::string &buf) {
    buf.clear();
    int c;
    while ((c = fpInTree.get()) != EOF) {
        if (c == '(' || c == ')' || c == ':' || c == ';' || c == ',') {
            /* standalone token */
            if (buf.empty() == 0) {
                buf += c;
                break;
            } else {
                fpInTree.unget();
                break;
            }
        } else if (std::isspace(c)) {
            if (!buf.empty()) {
                break;
            }
            /* else ignore whitespace at beginning of token */
        } else {
            /* not whitespace or standalone token */
            buf += c;
        }
    }
    return !buf.empty();
}

AbsNeighbourJoining(size_t)::TraversePostorder(int64_t node, std::vector<bool> &traversal, bool *pUp) {
    if (pUp != nullptr) {
        *pUp = false;
    }
    while (true) {
        assert(node >= 0);

        /* move to a child if possible */
        bool found = false;
        int iChild;
        for (iChild = 0; iChild < child[node].nChild; iChild++) {
            size_t childnode = child[node].child[iChild];
            if (!traversal[childnode]) {
                node = childnode;
                found = true;
                break;
            }
        }
        if (found) {
            continue; /* keep moving down */
        }
        if (!traversal[node]) {
            traversal[node] = true;
            return (node);
        }
        /* If we've already done this node, need to move up */
        if (node == root) {
            return -1; /* nowhere to go -- done traversing */
        }
        node = parent[node];
        /* If we go up to someplace that was already marked as visited, this is due
           to a change in topology, so return it marked as "up" */
        if (pUp != nullptr && traversal[node]) {
            *pUp = true;
            return node;
        }
    }
}

AbsNeighbourJoining(void)::readTreeError(const std::string &err, const std::string &token) {
    throw std::invalid_argument(strformat("Tree parse error: unexpected token '%s' -- %s",
                                          token.empty() ? "(End of file)" : token,
                                          err
    ));
}

#endif