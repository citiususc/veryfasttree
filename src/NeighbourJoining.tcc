
#ifndef FASTTREE_NEIGHBOURJOINING_TCC
#define FASTTREE_NEIGHBOURJOINING_TCC

#include "NeighbourJoining.h"
#include <list>
#include <unordered_map>
#include <sys/mman.h>
#include <cstring>
#include <fcntl.h>

#define AbsNeighbourJoining(...) \
template<typename Precision, template<class> class Operations> \
__VA_ARGS__ veryfasttree::NeighbourJoining<Precision, Operations>


AbsNeighbourJoining()::Profile::Profile(int64_t nPos, int64_t nConstraints) {
    typename op_t::Allocator alloc;
    diskLevel = 0;
    weights = alloc.allocate(nPos);
    codes = new char[nPos];
    vectors = nullptr;
    codeDist = nullptr;
    if (nConstraints == 0) {
        nOn = nullptr;
        nOff = nullptr;
    } else {
        nOn = new int64_t[nPos];
        nOff = new int64_t[nPos];
    }
    reset();
}

AbsNeighbourJoining()::Profile::Profile(int64_t nPos, int64_t nConstraints, uintptr_t &men, int nCodes, bool optAttr,
                                        bool test) {
    diskLevel = optAttr ? 2 : 1;
    size_t space = 10e6; // Just some big number.
    void *ptr = (void *) men;

    weights = (numeric_t *) std::align(op_t::ALIGNMENT, sizeof(Precision), ptr, space);
    ptr = weights + sizeof(Precision) * nPos;

    if (optAttr) {
        int nCodeSize = alignsz(nCodes, op_t::ALIGNMENT / sizeof(numeric_t));
        vectors = (numeric_t *) std::align(op_t::ALIGNMENT, sizeof(Precision), ptr, space);
        ptr = vectors + sizeof(Precision) * nPos * nCodeSize;

        codeDist = (numeric_t *) std::align(op_t::ALIGNMENT, sizeof(Precision), ptr, space);
        ptr = codeDist + sizeof(Precision) * nPos * nCodes;
    } else {
        vectors = nullptr;
        codeDist = nullptr;
    }

    if (nConstraints == 0) {
        nOn = nullptr;
        nOff = nullptr;
    } else {
        nOn = (int64_t *) ptr;
        ptr = nOn + nPos * sizeof(int64_t);

        nOff = (int64_t *) ptr;
        ptr = nOff + nPos * sizeof(int64_t);
    }

    codes = (char *) ptr;
    ptr = codes + nPos * sizeof(char);

    men = (uintptr_t) ptr;
    if (!test) {
        reset();
    }
}

AbsNeighbourJoining()::Profile::~Profile() {
    if (diskLevel < 1) {
        typename op_t::Allocator alloc;
        alloc.deallocate(weights, -1/*ignored*/);
        delete[] codes;
        if (nOn != nullptr) {
            delete[] nOn;
            delete[] nOff;
        }
    }
    reset();
}

AbsNeighbourJoining(void)::Profile::setVectorSize(size_t n, numeric_t val) {
    vectorsSize = n;
    if (diskLevel < 2) {
        typename op_t::Allocator alloc;
        if (vectors != nullptr) {
            alloc.deallocate(vectors, vectorsSize);
            vectors = nullptr;
        }
        if (n > 0) {
            vectors = alloc.allocate(n);
        }
    }
    for (size_t i = 0; i < n; i++) { vectors[i] = val; }
}

AbsNeighbourJoining(void)::Profile::resizeVector(size_t n) {
    if (vectors != nullptr && vectorsSize != n) {
        if (diskLevel < 2) {
            typename op_t::Allocator alloc;
            numeric_t *tmp = vectors;
            vectors = alloc.allocate(n);
            std::memcpy(vectors, tmp, n * sizeof(numeric_t));
            alloc.deallocate(tmp, vectorsSize);
        }
        vectorsSize = n;
    }
}


AbsNeighbourJoining(void)::Profile::setCodeDistSize(size_t n) {
    codeDistSize = n;
    if (diskLevel < 2) {
        typename op_t::Allocator alloc;
        if (codeDist != nullptr) {
            alloc.deallocate(codeDist, codeDistSize);
            codeDist = nullptr;
        }
        if (n > 0) {
            codeDist = alloc.allocate(n);
        }
    }
}

AbsNeighbourJoining(void)::Profile::reset() {
    nGaps = 0;
    nVectors = 0;
    setVectorSize(0, 0);
    setCodeDistSize(0);
}

AbsNeighbourJoining()::Rates::Rates(int64_t nRateCategories, int64_t nPos) {
    assert(nRateCategories >= 0);
    if (nRateCategories > 0) {
        rates.resize(nRateCategories, 1.0);
        ratecat.resize(nPos, 0);
    }
}

AbsNeighbourJoining(void)::Rates::reset(int64_t nRateCategories, int64_t nPos) {
    assert(nRateCategories >= 0);
    rates.clear();
    ratecat.clear();
    rates.resize(nRateCategories, 1.0);
    ratecat.resize(nPos, 0);
}

AbsNeighbourJoining()::TopHits::TopHits() {}

AbsNeighbourJoining()::TopHits::TopHits(const Options &options, int64_t _maxnodes, int64_t _m) :
        m(_m), q((int64_t) (0.5 + options.tophits2Mult * std::sqrt(m))),
        maxnodes(_maxnodes), topvisibleAge(0) {
    assert(m > 0);
    if (!options.useTopHits2nd || q >= m) {
        q = 0;
    }
    topHitsLists.resize(maxnodes, {{}, -1, 0});
    visible.resize(maxnodes, {-1, (numeric_t) 1e20});
    int64_t nTopVisible = (int64_t) (0.5 + options.topvisibleMult * m);
    topvisible.resize(nTopVisible, -1);
}

AbsNeighbourJoining()::
NeighbourJoining(Options &options, std::ostream &log, ProgressReport &progressReport,
                 std::vector<std::string> &seqs, int64_t nPos,
                 std::vector<std::string> &constraintSeqs,
                 DistanceMatrix <Precision, op_t::ALIGNMENT> &distanceMatrix,
                 TransitionMatrix <Precision, op_t::ALIGNMENT> &transmat) : log(log),
                                                                            options(options),
                                                                            progressReport(progressReport),
                                                                            nSeqs((int64_t) seqs.size()),
                                                                            nCodeSize(alignsz(options.nCodes,
                                                                                              op_t::ALIGNMENT /
                                                                                              sizeof(numeric_t))),
                                                                            distanceMatrix(distanceMatrix),
                                                                            transmat(transmat),
                                                                            constraintSeqs(constraintSeqs),
                                                                            outprofile(nPos, constraintSeqs.size()),
                                                                            rates(1, nPos) {
    this->root = -1;
    this->maxnode = nSeqs;
    this->nPos = nPos;
    this->maxnodes = 2 * nSeqs;


    seqsToProfiles(seqs);
    /* profiles from nSeq to maxnodes not yet exists */

    outProfile(outprofile, profiles, nSeqs);

    if (options.verbose > 10) {
        log << "Made out-profile" << std::endl;
    }

    totdiam = 0.0;

    diameter.resize(maxnodes, 0);
    varDiameter.resize(maxnodes, 0);
    selfdist.resize(maxnodes, 0);
    selfweight.resize(maxnodes);
    for (int64_t i = 0; i < (int64_t) nSeqs; i++) {
        selfweight[i] = (numeric_t)(nPos - profiles[i].nGaps);
        assert(selfweight[i] == (nPos - nGaps(i)));
    }

    outDistances.resize(maxnodes);
    nOutDistActive.resize(maxnodes, nSeqs * 10); /* unreasonably high value */
    // parent.empty()        /* so SetOutDistance ignores it */
    #pragma omp parallel for schedule(dynamic)
    for (int64_t i = 0; i < (int64_t) nSeqs; i++) {
        setOutDistance(i, nSeqs);
    }

    if (options.verbose > 2) {
        for (int64_t i = 0; i < 4 && i < (int64_t) nSeqs; i++) {
            log << strformat("Node %ld outdist %f", i, outDistances[i]) << std::endl;
        }
    }

    parent.resize(maxnodes, -1);
    branchlength.resize(maxnodes, 0); /* distance to parent */
    support.resize(maxnodes, -1.0);
    child.resize(maxnodes);
}

AbsNeighbourJoining()::~NeighbourJoining() {
    if (profilesMapSize > 0) {
        if (profilesMap != 0) {
            munmap((void *) profilesMap, profilesMapSize);
        }
        if (profilesFile != -1) {
            close(profilesFile);
        }
        remove(options.diskProfilesFile.c_str());
    }
}

AbsNeighbourJoining(void)::printDistances(std::vector<std::string> &names, std::ostream &out) {
    for (int64_t i = 0; i < (int64_t) nSeqs; i++) {
        std::cout << names[i];
        for (int64_t j = 0; j < (int64_t) nSeqs; j++) {
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

AbsNeighbourJoining(double)::totalLen() {
    double total_len = 0;
    for (int64_t iNode = 0; iNode < maxnode; iNode++) {
        total_len += std::fabs(branchlength[iNode]);
    }
    return total_len;
}

AbsNeighbourJoining(void)::branchlengthScale() {
    std::vector<numeric_t, typename op_t::Allocator> rates;
    std::vector<double> site_loglk;
    MLSiteRates(rates);
    MLSiteLikelihoodsByRate(rates, site_loglk);
    double scale = rescaleGammaLogLk(rates, site_loglk);

    for (int64_t i = 0; i < maxnodes; i++) {
        branchlength[i] *= scale;
    }
}

AbsNeighbourJoining(const std::vector<Precision, typename Operations<Precision>::Allocator> &)::getBranchlength() {
    return branchlength;
}

AbsNeighbourJoining(int64_t)::getMaxnode() {
    return maxnode;
}

AbsNeighbourJoining(int64_t)::getRateCategories() {
    return (int64_t) rates.rates.size();
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

/* Print topology using node indices as node names */
AbsNeighbourJoining(void)::printNJInternal(std::ostream &out, bool useLen) {
    if (nSeqs < 4) {
        return;
    }
    std::vector<std::pair<int64_t, int64_t>> stack(maxnodes);
    int64_t stackSize = 1;
    stack[0].first = root;
    stack[0].second = 0;

    while (stackSize > 0) {
        auto &last = stack[stackSize - 1];
        stackSize--;
        /* Save last, as we are about to overwrite it */
        int64_t node = last.first;
        int64_t end = last.second;

        if (node < (int64_t) nSeqs) {
            if (child[parent[node]].child[0] != node) {
                out << ",";
            }
            out << node;
            if (useLen) {
                out << strformat(":%.4f", branchlength[node]);
            }
        } else if (end) {
            out << strformat(")%ld", node);
            if (useLen) {
                out << strformat(":%.4f", branchlength[node]);
            }
        } else {
            if (node != root && child[parent[node]].child[0] != node) {
                out << ",";
            }
            out << "(";
            stackSize++;
            stack[stackSize - 1].first = node;
            stack[stackSize - 1].second = 1;
            Children &c = child[node];
            /* put children on in reverse order because we use the last one first */
            for (int i = c.nChild - 1; i >= 0; i--) {
                stackSize++;
                stack[stackSize - 1].first = c.child[i];
                stack[stackSize - 1].second = 0;
            }
        }
    }
    out << ";" << std::endl;
}

AbsNeighbourJoining(void)::seqsToProfiles(std::vector<std::string> &seqs) {
    profiles.reserve(maxnodes);
    profilesMapSize = 0;
    int64_t diskProfiles = (int64_t) std::ceil(maxnodes * options.diskProfilesRatio);
    if (diskProfiles > 0) {
        profilesMap = (uintptr_t) &diskProfiles;
        Profile(nPos, constraintSeqs.size(), profilesMap, options.nCodes, options.diskProfilesOpt,
                true); //Check profile size
        int64_t profileSize = (int64_t) (profilesMap - (uintptr_t) &diskProfiles) + op_t::ALIGNMENT;
        profilesMapSize = profileSize * diskProfiles;
        profilesMap = 0;

        profilesFile = open(options.diskProfilesFile.c_str(), O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
        if (profilesFile == -1) {
            throw std::runtime_error("disk memory file is invalid");
        }
        lseek(profilesFile, profilesMapSize, SEEK_SET);
        if (write(profilesFile, "", 1) == -1) {
            throw std::runtime_error("disk memory file truncation error");
        }
        void *men = mmap(0, profilesMapSize, PROT_READ | PROT_WRITE, MAP_PRIVATE, profilesFile, 0);
        if (men == MAP_FAILED) {
            throw std::runtime_error("profiles memory mapping fails");
        }
        profilesMap = (uintptr_t) men;

        for (int64_t i = 0; i < diskProfiles; i++) {
            profiles.emplace_back(nPos, constraintSeqs.size(), profilesMap, options.nCodes, options.diskProfilesOpt);
        }
    }
    for (int64_t i = (int64_t) profiles.size(); i < maxnodes; i++) {
        profiles.emplace_back(nPos, constraintSeqs.size());
    }
    uint8_t charToCode[256];
    std::unique_ptr<int64_t[]> counts(new int64_t[256 * options.threads]); /*Array of zeros*/

    for (int c = 0; c < 256; c++) {
        charToCode[c] = static_cast<uint8_t>(options.nCodes);
    }
    for (int i = 0; options.codesString[i]; i++) {
        charToCode[static_cast<int>(options.codesString[i])] = static_cast<uint8_t>(i);
    }
    charToCode['-'] = NOCODE;

    #pragma omp parallel
    {
        int64_t *lcounts = counts.get() + (256 * omp_get_thread_num());
        for (int c = 0; c < 256; c++) {
            lcounts[c] = 0;
        }
        #pragma omp for schedule(static)
        for (int64_t i = 0; i < (int64_t) nSeqs; i++) {
            auto &seq = seqs[i];
            auto &profile = profiles[i];
            for (int64_t j = 0; j < nPos; j++) {
                auto character = (uint8_t) seq[j];
                lcounts[character]++;
                auto c = charToCode[character];
                if (options.verbose > 10 && j < 2) {
                    log << strformat("pos %ld char %c code %u", j, seq[j], c) << std::endl;
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
            strrelease(seq);
            if (!constraintSeqs.empty()) {
                auto &constraintSeq = constraintSeqs[i];
                for (int64_t j = 0; j < (int64_t) constraintSeq.size(); j++) {
                    if (constraintSeq[j] == '1') {
                        profile.nOn[j] = 1;
                    } else if (constraintSeq[j] == '0') {
                        profile.nOff[j] = 1;
                    } else if (constraintSeq[j] != '-') {
                        #pragma omp critical
                        {
                            log << strformat("Constraint characters in unique sequence %ld replaced with gap: %c%ld",
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
        if (options.threads > 1) {
            #pragma omp for schedule(static)
            for (int i = 0; i < 256; i++) {
                for (int j = 1; j < options.threads; j++) {
                    counts[i] += counts[i + j * 256];
                }
            }
        }
    }
    int64_t totCount = 0;
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

AbsNeighbourJoining(int64_t)::activeAncestor(int64_t iNode) {
    if (iNode < 0) {
        return iNode;
    }
    while (parent[iNode] >= 0) {
        iNode = parent[iNode];
    }
    return iNode;
}

AbsNeighbourJoining(bool)::getVisible(int64_t nActive, TopHits &tophits, int64_t iNode, Besthit &visible) {
    if (iNode < 0 || parent[iNode] >= 0) {
        return false;
    }
    Hit &v = tophits.visible[iNode];
    if (v.j < 0 || parent[v.j] >= 0) {
        return false;
    }
    hitToBestHit(iNode, v, visible);
    setCriterion(nActive, visible);
    return true;
}

AbsNeighbourJoining(int64_t)::joinConstraintPenalty(int64_t node1, int64_t node2) {
    if (constraintSeqs.size() == 0) {
        return 0;
    }
    int64_t penalty = 0;
    for (int64_t iC = 0; iC < (int64_t) constraintSeqs.size(); iC++) {
        penalty += joinConstraintPenaltyPiece(node1, node2, iC);
    }
    return penalty;
}

AbsNeighbourJoining(int64_t)::joinConstraintPenaltyPiece(int64_t node1, int64_t node2, int64_t iC) {
    Profile &pOut = outprofile;
    Profile &p1 = profiles[node1];
    Profile &p2 = profiles[node2];
    int64_t nOn1 = p1.nOn[iC];
    int64_t nOff1 = p1.nOff[iC];
    int64_t nOn2 = p2.nOn[iC];
    int64_t nOff2 = p2.nOff[iC];
    int64_t nOnOut = pOut.nOn[iC] - nOn1 - nOn2;
    int64_t nOffOut = pOut.nOff[iC] - nOff1 - nOff2;

    if ((nOn1 + nOff1) > 0 && (nOn2 + nOff2) > 0 && (nOnOut + nOffOut) > 0) {
        /* code is -1 for split, 0 for off, 1 for on */
        int64_t code1 = (nOn1 > 0 && nOff1 > 0) ? -1 : (nOn1 > 0 ? 1 : 0);
        int64_t code2 = (nOn2 > 0 && nOff2 > 0) ? -1 : (nOn2 > 0 ? 1 : 0);
        int64_t code3 = (nOnOut > 0 && nOffOut > 0) ? -1 : (nOnOut > 0 ? 1 : 0);
        int64_t nSplit = (code1 == -1 ? 1 : 0) + (code2 == -1 ? 1 : 0) + (code3 == -1 ? 1 : 0);
        int64_t nOn = (code1 == 1 ? 1 : 0) + (code2 == 1 ? 1 : 0) + (code3 == 1 ? 1 : 0);
        if (nSplit == 1 && nOn == 1) {
            return splitConstraintPenalty(nOn1 + nOn2, nOff1 + nOff2, nOnOut, nOffOut);
        }
    }
    /* else */
    return 0;
}

/* Minimum number of constrained leaves that need to be moved
   to satisfy the constraint (or 0 if constraint is satisfied)
   Defining it this way should ensure that SPR moves that break
   constraints get a penalty
*/
AbsNeighbourJoining(int64_t)::splitConstraintPenalty(int64_t nOn1, int64_t nOff1, int64_t nOn2, int64_t nOff2) {
    return (nOn1 + nOff2 < nOn2 + nOff1 ?
            (nOn1 < nOff2 ? nOn1 : nOff2)
                                        : (nOn2 < nOff1 ? nOn2 : nOff1));
}

/* Computes support for (A,B),(C,D) compared to that for (A,C),(B,D) and (A,D),(B,C) */
AbsNeighbourJoining(double)::
splitSupport(Profile &pA, Profile &pB, Profile &pC, Profile &pD, const std::vector<int64_t> &col) {
    /* Note distpieces are weighted */
    std::vector<double> distpieces(6 * nPos);
    std::vector<double> weights(6 * nPos);

    int64_t iFreqA = 0;
    int64_t iFreqB = 0;
    int64_t iFreqC = 0;
    int64_t iFreqD = 0;
    for (int64_t i = 0; i < nPos; i++) {
        numeric_t *fA = getFreq(pA, i, iFreqA);
        numeric_t *fB = getFreq(pB, i, iFreqB);
        numeric_t *fC = getFreq(pC, i, iFreqC);
        numeric_t *fD = getFreq(pD, i, iFreqD);

        weights[qAB * nPos + i] = pA.weights[i] * pB.weights[i];
        weights[qAC * nPos + i] = pA.weights[i] * pC.weights[i];
        weights[qAD * nPos + i] = pA.weights[i] * pD.weights[i];
        weights[qBC * nPos + i] = pB.weights[i] * pC.weights[i];
        weights[qBD * nPos + i] = pB.weights[i] * pD.weights[i];
        weights[qCD * nPos + i] = pC.weights[i] * pD.weights[i];

        distpieces[qAB * nPos + i] =
                weights[qAB * nPos + i] * profileDistPiece(pA.codes[i], pB.codes[i], fA, fB, nullptr);
        distpieces[qAC * nPos + i] =
                weights[qAC * nPos + i] * profileDistPiece(pA.codes[i], pC.codes[i], fA, fC, nullptr);
        distpieces[qAD * nPos + i] =
                weights[qAD * nPos + i] * profileDistPiece(pA.codes[i], pD.codes[i], fA, fD, nullptr);
        distpieces[qBC * nPos + i] =
                weights[qBC * nPos + i] * profileDistPiece(pB.codes[i], pC.codes[i], fB, fC, nullptr);
        distpieces[qBD * nPos + i] =
                weights[qBD * nPos + i] * profileDistPiece(pB.codes[i], pD.codes[i], fB, fD, nullptr);
        distpieces[qCD * nPos + i] =
                weights[qCD * nPos + i] * profileDistPiece(pC.codes[i], pD.codes[i], fC, fD, nullptr);
    }
    assert(iFreqA == pA.nVectors);
    assert(iFreqB == pB.nVectors);
    assert(iFreqC == pC.nVectors);
    assert(iFreqD == pD.nVectors);

    double totpieces[6];
    double totweights[6];
    double dists[6];
    for (int j = 0; j < 6; j++) {
        totpieces[j] = 0.0;
        totweights[j] = 0.0;
        for (int64_t i = 0; i < nPos; i++) {
            totpieces[j] += distpieces[j * nPos + i];
            totweights[j] += weights[j * nPos + i];
        }
        dists[j] = totweights[j] > 0.01 ? totpieces[j] / totweights[j] : 3.0;
        if (options.logdist) {
            dists[j] = logCorrect(dists[j]);
        }
    }

    /* Support1 = Support(AB|CD over AC|BD) = d(A,C)+d(B,D)-d(A,B)-d(C,D)
       Support2 = Support(AB|CD over AD|BC) = d(A,D)+d(B,C)-d(A,B)-d(C,D)
    */
    double support1 = dists[qAC] + dists[qBD] - dists[qAB] - dists[qCD];
    double support2 = dists[qAD] + dists[qBC] - dists[qAB] - dists[qCD];

    if (support1 < 0 || support2 < 0) {
        options.debug.nSuboptimalSplits++;    /* Another split seems superior */
    }

    assert(options.nBootstrap > 0);
    int64_t nSupport = 0;

    for (int64_t iBoot = 0; iBoot < options.nBootstrap; iBoot++) {
        const int64_t *colw = &col[nPos * iBoot];

        for (int j = 0; j < 6; j++) {
            double totp = 0;
            double totw = 0;
            double *d = &distpieces[j * nPos];
            double *w = &weights[j * nPos];
            for (int64_t i = 0; i < nPos; i++) {
                int64_t c = colw[i];
                totp += d[c];
                totw += w[c];
            }
            dists[j] = totw > 0.01 ? totp / totw : 3.0;
            if (options.logdist) {
                dists[j] = logCorrect(dists[j]);
            }
        }
        support1 = dists[qAC] + dists[qBD] - dists[qAB] - dists[qCD];
        support2 = dists[qAD] + dists[qBC] - dists[qAB] - dists[qCD];
        if (support1 > 0 && support2 > 0) {
            nSupport++;
        }
    } /* end loop over bootstrap replicates */
    return nSupport / (double) options.nBootstrap;
}

AbsNeighbourJoining(void)::resampleColumns(std::vector<int64_t> &col) {
    col.resize(nPos * options.nBootstrap);
    for (int64_t i = 0; i < options.nBootstrap; i++) {
        for (int64_t j = 0; j < nPos; j++) {
            int64_t pos = (int64_t) (knuth_rand() * nPos);
            if (pos < 0) {
                pos = 0;
            } else if (pos == nPos) {
                pos = nPos - 1;
            }
            col[i * nPos + j] = pos;
        }
    }
    if (options.verbose > 5) {
        for (int i = 0; i < 3 && i < options.nBootstrap; i++) {
            log << "Boot" << i;
            for (int64_t j = 0; j < nPos; j++) {
                log << "\t" << col[i * nPos + j];
            }
            log << std::endl;
        }
    }
}

AbsNeighbourJoining(template<typename Profile_t> void)::
outProfile(Profile &out, std::vector<Profile_t> &_profiles, int64_t nProfiles) {
    out.reset();
    double inweight = 1.0 / (double) nProfiles;   /* The maximal output weight is 1.0 */

    #pragma omp parallel
    {
        /* First, set weights -- code is always NOCODE, prevent weight=0 */
        #pragma omp for schedule(dynamic)
        for (int64_t i = 0; i < nPos; i++) {
            out.weights[i] = 0;
            for (int64_t in = 0; in < nProfiles; in++) {
                out.weights[i] += asRef(_profiles[in]).weights[i] * inweight;
            }
            if (out.weights[i] <= 0) {
                out.weights[i] = (numeric_t) 1e-20;
            } /* always store a vector */
            out.codes[i] = NOCODE;        /* outprofile is normally complicated */
        }

        #pragma omp single
        {
            out.nVectors = nPos;
            /* Initialize the frequencies to 0 */
            out.setVectorSize(out.nVectors * nCodeSize, 0);
        }
        #pragma omp barrier

        std::vector<numeric_t, typename op_t::Allocator> tmp(out.nVectors * nCodeSize, 0);
        #pragma omp for schedule(dynamic)
        for (int64_t in = 0; in < nProfiles; in++) {
            int64_t iFreqOut = 0;
            int64_t iFreqIn = 0;
            for (int64_t i = 0; i < nPos; i++) {
                numeric_t *fIn = getFreq(asRef(_profiles[in]), i, iFreqIn);
                getFreq(out, i, iFreqOut);
                if (asRef(_profiles[in]).weights[i] > 0) {
                    numeric_t *lfOut = &tmp[nCodeSize * (iFreqOut - 1)];
                    addToFreq(lfOut, asRef(_profiles[in]).weights[i], asRef(_profiles[in]).codes[i], fIn);
                }
            }
            assert(iFreqOut == out.nVectors);
            assert(iFreqIn == asRef(_profiles[in]).nVectors);
        }

        #pragma omp critical
        {
            operations.vector_add(&out.vectors[0], &tmp[0], (int64_t) out.vectorsSize);
        }
        #pragma omp barrier

        /* And normalize the frequencies to sum to 1 */
        int64_t iFreqOut = 0;
        #pragma omp for schedule(dynamic)
        for (int64_t i = 0; i < nPos; i++) {
            iFreqOut = i;// code is always NOCODE and weight>0
            numeric_t *fOut = getFreq(out, i, iFreqOut);
            assert(fOut != nullptr);
            normalizeFreq(fOut, distanceMatrix);
        }

        assert(nPos == out.nVectors);
        if (options.verbose > 10) {
            log << strformat("Average %ld profiles", nProfiles) << std::endl;
        }

        if (distanceMatrix) {
            setCodeDist(out, /*shared*/true);
        }

        /* Compute constraints */
        #pragma omp for schedule(static)
        for (int64_t i = 0; i < (int64_t) constraintSeqs.size(); i++) {
            for (int64_t in = 0; in < nProfiles; in++) {
                out.nOn[i] += asRef(_profiles[in]).nOn[i];
                out.nOff[i] += asRef(_profiles[in]).nOff[i];
            }
        }
    }

}

AbsNeighbourJoining(void)::addToFreq(numeric_t fOut[], double weight, int64_t codeIn, numeric_t fIn[]) {
    addToFreq(fOut, weight, codeIn, fIn, distanceMatrix);
}

AbsNeighbourJoining(void)::addToFreq(numeric_t fOut[], double weight, int64_t codeIn, numeric_t fIn[],
                                     DistanceMatrix <Precision, op_t::ALIGNMENT> &dmat) {
    assert(fOut != nullptr);
    if (fIn != nullptr) {
        operations.vector_add_mult(fOut, fIn, weight, options.nCodes);
    } else if (dmat) {
        assert(codeIn != NOCODE);
        operations.vector_add_mult(fOut, dmat.codeFreq[codeIn], weight, options.nCodes);
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
    normalizeFreq(freq, distanceMatrix);
}

AbsNeighbourJoining(void)::normalizeFreq(numeric_t freq[], DistanceMatrix <Precision, op_t::ALIGNMENT> &dmat) {
    double total_freq = 0;
    if (dmat) {
        /* The total frequency is dot_product(true_frequencies, 1)
           So we rotate the 1 vector by eigeninv (stored in eigentot)
        */
        total_freq = operations.vector_multiply_sum(freq, dmat.eigentot, options.nCodes);
    } else {
        for (int k = 0; k < options.nCodes; k++) {
            total_freq += freq[k];
        }
    }
    if (total_freq > options.fPostTotalTolerance) {
        numeric_t inverse_weight = 1.0 / total_freq;
        operations.vector_multiply_by(freq, inverse_weight, options.nCodes, freq);
    } else {
        /* This can happen if we are in a very low-weight region, e.g. if a mostly-gap position gets weighted down
           repeatedly; just set them all to arbitrary but legal values */
        if (!dmat) {
            for (int k = 0; k < options.nCodes; k++) {
                freq[k] = 1.0 / options.nCodes;
            }
        } else {
            for (int k = 0; k < options.nCodes; k++) {
                freq[k] = dmat.codeFreq[0][k];
            }
        }
    }
}

AbsNeighbourJoining(void)::setCodeDist(Profile &profile, bool shared) {
    if (shared) {
        int64_t iFreq = 0;
        #pragma omp for schedule(dynamic)
        for (int64_t i = 0; i < nPos; i++) {
            numeric_t *f = getFreq(profile, i, iFreq);
            assert(f != nullptr);
            for (int k = 0; k < options.nCodes; k++) {
                profile.codeDist[i * options.nCodes + k] = profileDistPiece(profile.codes[i], k, f, nullptr, nullptr);
            }
        }
    } else {
        if (profile.codeDistSize == 0) {
            profile.setCodeDistSize(nPos * options.nCodes);
        }
        int64_t iFreq = 0;
        for (int64_t i = 0; i < nPos; i++) {
            numeric_t *f = getFreq(profile, i, iFreq);

            for (int k = 0; k < options.nCodes; k++) {
                profile.codeDist[i * options.nCodes + k] = profileDistPiece(profile.codes[i], k, f, nullptr, nullptr);
            }
        }
        assert(iFreq == profile.nVectors);
    }
}

AbsNeighbourJoining(double)::
profileDistPiece(int64_t code1, int64_t code2, numeric_t f1[], numeric_t f2[], numeric_t codeDist2[]) {
    if (distanceMatrix) {
        if (code1 != NOCODE && code2 != NOCODE) { /* code1 vs code2 */
            return distanceMatrix.distances[code1][code2];
        } else if (codeDist2 != nullptr && code1 != NOCODE) { /* code1 vs. codeDist2 */
            return codeDist2[code1];
        } else { /* f1 vs f2 */
            if (f1 == nullptr) {
                if (code1 == NOCODE) { return (10.0); }
                f1 = &distanceMatrix.codeFreq[code1][0];
            }
            if (f2 == nullptr) {
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
                if (f2 == nullptr) { return (10.0); }
                return 1.0 - f2[code1]; /* code1 vs. f2 */
            }
        } else {
            if (code2 != NOCODE) {
                if (f1 == nullptr) { return (10.0); }
                return 1.0 - f1[code2]; /* f1 vs code2 */
            } else { /* f1 vs. f2 */
                if (f1 == nullptr || f2 == nullptr) { return (10.0); }
                double piece = 1.0;
                for (int k = 0; k < options.nCodes; k++) {
                    piece -= f1[k] * f2[k];
                }
                return piece;
            }
        }
    }
}

AbsNeighbourJoining(void)::updateOutProfile(Profile &out, Profile &old1, Profile &old2, Profile &_new,
                                            int64_t nActiveOld) {
    int64_t iFreqOut = 0;
    int64_t iFreq1 = 0;
    int64_t iFreq2 = 0;
    int64_t iFreqNew = 0;
    assert(nActiveOld > 0);

    for (int64_t i = 0; i < nPos; i++) {
        numeric_t *fOut = getFreq(out, i, iFreqOut);
        numeric_t *fOld1 = getFreq(old1, i, iFreq1);
        numeric_t *fOld2 = getFreq(old2, i, iFreq2);
        numeric_t *fNew = getFreq(_new, i, iFreqNew);

        assert(out.codes[i] == NOCODE && fOut != nullptr); /* No no-vector optimization for outprofiles */
        if (options.verbose > 3 && i < 3) {
            log << strformat("Updating out-profile position %ld weight %f (mult %f)",
                             i, out.weights[i], out.weights[i] * nActiveOld) << std::endl;
        }
        double originalMult = out.weights[i] * nActiveOld;
        double newMult = originalMult + _new.weights[i] - old1.weights[i] - old2.weights[i];
        out.weights[i] = newMult / (nActiveOld - 1);
        if (out.weights[i] <= 0) {
            out.weights[i] = (numeric_t) 1e-20; /* always use the vector */
        }

        for (int64_t k = 0; k < options.nCodes; k++) {
            fOut[k] *= originalMult;
        }

        if (old1.weights[i] > 0) {
            addToFreq(fOut, -old1.weights[i], old1.codes[i], fOld1);
        }
        if (old2.weights[i] > 0) {
            addToFreq(fOut, -old2.weights[i], old2.codes[i], fOld2);
        }
        if (_new.weights[i] > 0) {
            addToFreq(fOut, _new.weights[i], _new.codes[i], fNew);
        }

        /* And renormalize */
        normalizeFreq(fOut);

        if (options.verbose > 2 && i < 3) {
            log << strformat("Updated out-profile position %ld weight %f (mult %f)",
                             i, out.weights[i], out.weights[i] * nActiveOld);
            if (out.weights[i] > 0) {
                for (int k = 0; k < options.nCodes; k++) {
                    log << strformat(" %c:%f", (distanceMatrix ? '?' : options.codesString[k]), fOut[k]);
                }
            }
            log << std::endl;
        }
    }
    assert(iFreqOut == out.nVectors);
    assert(iFreq1 == old1.nVectors);
    assert(iFreq2 == old2.nVectors);
    assert(iFreqNew == _new.nVectors);
    if (distanceMatrix) {
        setCodeDist(out);
    }

    /* update constraints -- note in practice this should be a no-op */
    for (int64_t i = 0; i < (int64_t) constraintSeqs.size(); i++) {
        out.nOn[i] += _new.nOn[i] - old1.nOn[i] - old2.nOn[i];
        out.nOff[i] += _new.nOff[i] - old1.nOff[i] - old2.nOff[i];
    }
}

AbsNeighbourJoining(void)::setOutDistance(int64_t iNode, int64_t nActive) {
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
            log << strformat("NewOutDist for %ld %f from dist %f selfd %f diam %f totdiam %f newActive %ld",
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
            for (int64_t j = 0; j < maxnode; j++) {
                if (j != iNode && (parent.empty() || parent[j] < 0)) {
                    Besthit bh;
                    profileDist(profiles[iNode], profiles[j], bh);
                    total_pd += bh.dist;
                    total += bh.dist - (diameter[iNode] + diameter[j]);
                }
            }
            log << strformat("OutDist for Node %ld %f truth %f profiled %f truth %f pd_err %f",
                             iNode, outDistances[iNode], total, pdistOutWithoutA, total_pd,
                             fabs(pdistOutWithoutA - total_pd)) << std::endl;
        };
    }
}

AbsNeighbourJoining(void)::setCriterion(int64_t nActive, Besthit &join) {
    if (join.i < 0 || join.j < 0 || parent[join.i] >= 0 || parent[join.j] >= 0) {
        return;
    }
    assert(nOutDistActive[join.i] >= nActive);
    assert(nOutDistActive[join.j] >= nActive);

    int64_t nDiffAllow = options.tophitsMult > 0 ? (int64_t) (nActive * options.staleOutLimit) : 0;
    if (nOutDistActive[join.i] - nActive > nDiffAllow) {
        setOutDistance(join.i, nActive);
    }
    if (nOutDistActive[join.j] - nActive > nDiffAllow) {
        setOutDistance(join.j, nActive);
    }
    double outI = outDistances[join.i];
    if (nOutDistActive[join.i] != nActive) {
        outI *= (nActive - 1) / (double) (nOutDistActive[join.i] - 1);
    }
    double outJ = outDistances[join.j];
    if (nOutDistActive[join.j] != nActive) {
        outJ *= (nActive - 1) / (double) (nOutDistActive[join.j] - 1);
    }
    join.criterion = join.dist - (outI + outJ) / (double) (nActive - 2);
    if (options.verbose > 2 && nActive <= 5) {
        log << strformat("Set Criterion to join %ld %ld with nActive=%ld dist+penalty %.3f criterion %.3f",
                         join.i, join.j, nActive, join.dist, join.criterion) << std::endl;
    }
}

AbsNeighbourJoining(void)::setDistCriterion(int64_t nActive, Besthit &hit) {
    if (hit.i < (int64_t) nSeqs && hit.j < (int64_t) nSeqs) {
        seqDist(profiles[hit.i].codes, profiles[hit.j].codes, hit);
    } else {
        profileDist(profiles[hit.i], profiles[hit.j], hit);
        hit.dist -= (diameter[hit.i] + diameter[hit.j]);
    }
    hit.dist += options.constraintWeight * (double) joinConstraintPenalty(hit.i, hit.j);
    setCriterion(nActive, hit);
}

AbsNeighbourJoining(double)::SHSupport(const std::vector<int64_t> &col, double loglk[3],
                                       std::vector<double> &site_likelihoods) {
    assert(options.nBootstrap > 0);
    double delta1 = loglk[0] - loglk[1];
    double delta2 = loglk[0] - loglk[2];
    double delta = delta1 < delta2 ? delta1 : delta2;

    std::vector<double> siteloglk(3 * nPos);
    for (int i = 0; i < 3 * nPos; i++) {
        siteloglk[i] = std::log(site_likelihoods[i]);
    }

    int64_t nSupport = 0;
    for (int64_t iBoot = 0; iBoot < options.nBootstrap; iBoot++) {
        double resampled[3];
        for (int i = 0; i < 3; i++) {
            resampled[i] = -loglk[i];
        }
        for (int64_t j = 0; j < nPos; j++) {
            int64_t pos = col[iBoot * nPos + j];
            for (int i = 0; i < 3; i++) {
                resampled[i] += siteloglk[i * nPos + pos];
            }
        }
        int iBest = 0;
        for (int i = 1; i < 3; i++) {
            if (resampled[i] > resampled[iBest]) {
                iBest = i;
            }
        }
        double resample1 = resampled[iBest] - resampled[(iBest + 1) % 3];
        double resample2 = resampled[iBest] - resampled[(iBest + 2) % 3];
        double resampleDelta = resample1 < resample2 ? resample1 : resample2;
        if (resampleDelta < delta) {
            nSupport++;
        }
    }

    return nSupport / (double) options.nBootstrap;
}

AbsNeighbourJoining(void)::profileDist(Profile &profile1, Profile &profile2, Besthit &hit) {
    double top = 0;
    double denom = 0;
    int64_t iFreq1 = 0;
    int64_t iFreq2 = 0;
    if (!omp_in_parallel() && options.threads > 1 && options.threadsLevel > 0) {
        std::vector<int64_t> hasFreq1(nPos, 0);
        std::vector<int64_t> hasFreq2(nPos, 0);
        std::vector<int64_t> found1(options.threads, 0);
        std::vector<int64_t> found2(options.threads, 0);

        #pragma omp parallel firstprivate(iFreq1, iFreq2)
        {
            double top2 = top;
            double denom2 = denom;
            int id = omp_get_thread_num();
            #pragma omp for schedule(static, (nPos / options.threads) + 1)
            for (int64_t i = 0; i < nPos; i++) {
                hasFreq1[i] = getFreq(profile1, i, iFreq1) != nullptr ? 1 : 0;
                hasFreq2[i] = getFreq(profile2, i, iFreq2) != nullptr ? 1 : 0;
            }

            found1[id] = iFreq1;
            found2[id] = iFreq2;
            #pragma omp barrier
            iFreq1 = 0;
            iFreq2 = 0;
            for (int i = 0; i < id; i++) {
                iFreq1 += found1[i];
                iFreq2 += found2[i];
            }

            #pragma omp for schedule(static, (nPos / options.threads) + 1)
            for (int64_t i = 0; i < nPos; i++) {
                if (hasFreq1[i] > 0) {
                    hasFreq1[i] = iFreq1++;
                }
                if (hasFreq2[i] > 0) {
                    hasFreq2[i] = iFreq2++;
                }
            }

            assert(found1.back() == profile1.nVectors);
            assert(found1.back() == profile2.nVectors);

            #pragma omp for schedule(dynamic)
            for (int64_t i = 0; i < nPos; i++) {
                if (profile1.weights[i] > 0 && profile2.weights[i] > 0) {
                    numeric_t *f1 = getFreq(profile1, i, hasFreq1[i]);
                    numeric_t *f2 = getFreq(profile2, i, hasFreq2[i]);
                    double weight = profile1.weights[i] * profile2.weights[i];
                    denom2 += weight;
                    double piece = profileDistPiece(profile1.codes[i], profile2.codes[i], f1, f2,
                                                    (!profile2.codeDistSize == 0 ? &profile2.codeDist[i *
                                                                                                      options.nCodes]
                                                                                 : nullptr));
                    top2 += weight * piece;
                }
            }

            #pragma omp critical
            {
                denom += denom2;
                top += top2;
            }
        }
    } else {
        for (int64_t i = 0; i < nPos; i++) {
            numeric_t *f1 = getFreq(profile1, i, iFreq1);
            numeric_t *f2 = getFreq(profile2, i, iFreq2);
            if (profile1.weights[i] > 0 && profile2.weights[i] > 0) {
                double weight = profile1.weights[i] * profile2.weights[i];
                denom += weight;
                double piece = profileDistPiece(profile1.codes[i], profile2.codes[i], f1, f2,
                                                (!profile2.codeDistSize == 0 ? &profile2.codeDist[i * options.nCodes]
                                                                             : nullptr));
                top += weight * piece;
            }
        }
        assert(iFreq1 == profile1.nVectors);
        assert(iFreq2 == profile2.nVectors);
    }
    hit.weight = denom > 0 ? denom : 0.01; /* 0.01 is an arbitrarily low value of weight (normally >>1) */
    hit.dist = denom > 0 ? top / denom : 1;
    options.debug.profileOps++;
}

AbsNeighbourJoining(double)::pairLogLk(Profile &p1, Profile &p2, double length, double site_likelihoods[]) {
    double lk = 1.0;
    double loglk = 0.0;        /* stores underflow of lk during the loop over positions */

    assert(!rates.rates.empty());
    std::vector<numeric_t, typename op_t::Allocator> expeigenRates;
    if (transmat) {
        expEigenRates(length, expeigenRates);
    }

    if (!transmat) {    /* Jukes-Cantor */
        assert (options.nCodes == 4);
        std::vector<double> pSame;
        std::vector<double> pDiff;

        pSameVector(length, pSame);
        pDiffVector(pSame, pDiff);

        int64_t iFreqA = 0;
        int64_t iFreqB = 0;
        for (int64_t i = 0; i < nPos; i++) {
            int64_t iRate = rates.ratecat[i];
            double wA = p1.weights[i];
            double wB = p2.weights[i];
            int64_t codeA = p1.codes[i];
            int64_t codeB = p2.codes[i];
            numeric_t *fA = getFreq(p1, i, iFreqA);
            numeric_t *fB = getFreq(p2, i, iFreqB);
            double lkAB = 0;

            if (fA == nullptr && fB == nullptr) {
                if (codeA == NOCODE) {    /* A is all gaps */
                    /* gap to gap is sum(j) 0.25 * (0.25 * pSame + 0.75 * pDiff) = sum(i) 0.25*0.25 = 0.25
                       gap to any character gives the same result
                    */
                    lkAB = 0.25;
                } else if (codeB == NOCODE) { /* B is all gaps */
                    lkAB = 0.25;
                } else if (codeA == codeB) { /* A and B match */
                    lkAB = pSame[iRate] * wA * wB + 0.25 * (1 - wA * wB);
                } else {        /* codeA != codeB */
                    lkAB = pDiff[iRate] * wA * wB + 0.25 * (1 - wA * wB);
                }
            } else if (fA == nullptr) {
                /* Compare codeA to profile of B */
                if (codeA == NOCODE)
                    lkAB = 0.25;
                else
                    lkAB = wA * (pDiff[iRate] + fB[codeA] * (pSame[iRate] - pDiff[iRate])) + (1.0 - wA) * 0.25;
                /* because lkAB = wA * P(codeA->B) + (1-wA) * 0.25
                   P(codeA -> B) = sum(j) P(B==j) * (j==codeA ? pSame : pDiff)
                   = sum(j) P(B==j) * pDiff +
                   = pDiff + P(B==codeA) * (pSame-pDiff)
                */
            } else if (fB == nullptr) { /* Compare codeB to profile of A */
                if (codeB == NOCODE) {
                    lkAB = 0.25;
                } else {
                    lkAB = wB * (pDiff[iRate] + fA[codeB] * (pSame[iRate] - pDiff[iRate])) + (1.0 - wB) * 0.25;
                }
            } else { /* both are full profiles */
                for (int j = 0; j < 4; j++) {
                    lkAB += fB[j] * (fA[j] * pSame[iRate] + (1 - fA[j]) * pDiff[iRate]); /* P(A|B) */
                }
            }
            assert(lkAB > 0);
            lk *= lkAB;
            while (lk < Constants::LkUnderflow) {
                lk *= Constants::LkUnderflowInv;
                loglk -= Constants::LogLkUnderflow;
            }
            if (site_likelihoods != nullptr) {
                site_likelihoods[i] *= lkAB;
            }
        }
    } else if (options.nCodes == 4) {    /* matrix model on nucleotides */
        int64_t iFreqA = 0;
        int64_t iFreqB = 0;
        numeric_t fAmix[4], fBmix[4];
        numeric_t *fGap = &transmat.codeFreq[NOCODE][0];

        for (int64_t i = 0; i < nPos; i++) {
            int64_t iRate = rates.ratecat[i];
            numeric_t *expeigen = &expeigenRates[iRate * nCodeSize];
            double wA = p1.weights[i];
            double wB = p2.weights[i];
            if (wA == 0 && wB == 0 && p1.codes[i] == NOCODE && p2.codes[i] == NOCODE) {
                /* Likelihood of A vs B is 1, so nothing changes
                   Do not need to advance iFreqA or iFreqB */
                continue;
            }
            numeric_t *fA = getFreq(p1, i, iFreqA);
            numeric_t *fB = getFreq(p2, i, iFreqB);
            if (fA == nullptr) {
                fA = &transmat.codeFreq[(int) p1.codes[i]][0];
            }
            if (wA > 0.0 && wA < 1.0) {
                for (int64_t j = 0; j < 4; j++) {
                    fAmix[j] = wA * fA[j] + (1.0 - wA) * fGap[j];
                }
                fA = fAmix;
            }
            if (fB == nullptr) {
                fB = &transmat.codeFreq[(int) p2.codes[i]][0];
            }
            if (wB > 0.0 && wB < 1.0) {
                for (int64_t j = 0; j < 4; j++) {
                    fBmix[j] = wB * fB[j] + (1.0 - wB) * fGap[j];
                }
                fB = fBmix;
            }
            /* SSE3 instructions do not speed this step up:
           numeric_t lkAB = vector_multiply3_sum(expeigen, fA, fB); */
            // dsp this is where check for <=0 was added in 2.1.1.LG
            double lkAB = 0;
            for (int j = 0; j < 4; j++) {
                lkAB += expeigen[j] * fA[j] * fB[j];
            }
            assert(lkAB > 0);
            if (site_likelihoods != nullptr)
                site_likelihoods[i] *= lkAB;
            lk *= lkAB;
            while (lk < Constants::LkUnderflow) {
                lk *= Constants::LkUnderflowInv;
                loglk -= Constants::LogLkUnderflow;
            }
            while (lk > Constants::LkUnderflowInv) {
                lk *= Constants::LkUnderflow;
                loglk += Constants::LogLkUnderflow;
            }
        }
    } else if (options.nCodes == 20) {    /* matrix model on amino acids */
        int64_t iFreqA = 0;
        int64_t iFreqB = 0;
        numeric_t fAmix[20], fBmix[20];
        numeric_t *fGap = &transmat.codeFreq[NOCODE][0];

        for (int64_t i = 0; i < nPos; i++) {
            int64_t iRate = rates.ratecat[i];
            numeric_t *expeigen = &expeigenRates[iRate * nCodeSize];
            double wA = p1.weights[i];
            double wB = p2.weights[i];
            if (wA == 0 && wB == 0 && p1.codes[i] == NOCODE && p2.codes[i] == NOCODE) {
                /* Likelihood of A vs B is 1, so nothing changes
                   Do not need to advance iFreqA or iFreqB */
                continue;
            }
            numeric_t *fA = getFreq(p1, i, iFreqA);
            numeric_t *fB = getFreq(p2, i, iFreqB);
            if (fA == nullptr) {
                fA = &transmat.codeFreq[(int) p1.codes[i]][0];
            }
            if (wA > 0.0 && wA < 1.0) {
                for (int j = 0; j < 20; j++) {
                    fAmix[j] = wA * fA[j] + (1.0 - wA) * fGap[j];
                }
                fA = fAmix;
            }
            if (fB == nullptr) {
                fB = &transmat.codeFreq[(int) p2.codes[i]][0];
            }
            if (wB > 0.0 && wB < 1.0) {
                for (int j = 0; j < 20; j++) {
                    fBmix[j] = wB * fB[j] + (1.0 - wB) * fGap[j];
                }
                fB = fBmix;
            }
            numeric_t lkAB = operations.vector_multiply3_sum(expeigen, fA, fB, 20);
            if (!(lkAB > 0)) {
                /* If this happens, it indicates a numerical problem that needs to be addressed elsewhere,
                   so report all the details */
                log << "# VeryFastTree.c::PairLogLk -- numerical problem!" << std::endl;
                log << "# This block is intended for loading into R" << std::endl;

                log << strformat("lkAB = %.8g", lkAB) << std::endl;
                log << strformat(
                        "Branch_length= %.8g\nalignment_position=%ld\nnCodes=%d\nrate_category=%ld\nrate=%.8g",
                        length, i, options.nCodes, iRate, rates.rates[iRate]) << std::endl;
                log << strformat("wA=%.8g\nwB=%.8g", wA, wB) << std::endl;
                log << strformat("codeA = %d\ncodeB = %d", p1.codes[i], p2.codes[i]) << std::endl;

                log << "fA = c(";
                for (int j = 0; j < options.nCodes; j++) {
                    log << (j == 0 ? "" : ",") << strformat(" %.8g", fA[j]);
                }
                log << ")" << std::endl;

                log << "fB = c(";
                for (int j = 0; j < options.nCodes; j++) {
                    log << (j == 0 ? "" : ",") << strformat(" %.8g", fB[j]);
                }
                log << ")" << std::endl;

                log << "stat = c(";
                for (int j = 0; j < options.nCodes; j++) {
                    log << (j == 0 ? "" : ",") << strformat(" %.8g", transmat.stat[j]);
                }
                log << ")" << std::endl;

                log << "eigenval = c(";
                for (int j = 0; j < options.nCodes; j++) {
                    log << (j == 0 ? "" : ",") << strformat(" %.8g", transmat.eigenval[j]);
                }
                log << ")" << std::endl;

                log << "expeigen = c(";
                for (int j = 0; j < options.nCodes; j++) {
                    log << (j == 0 ? "" : ",") << strformat(" %.8g", expeigen[j]);
                }
                log << ")" << std::endl;

                log << "codeFreq = c(";
                for (int j = 0; j < options.nCodes; j++) {
                    for (int k = 0; k < options.nCodes; k++) {
                        log << (j == 0 && k == 0 ? "" : ",") << strformat(" %.8g", transmat.codeFreq[j][k]);
                    }
                }
                log << ")" << std::endl;

                log << "eigeninv = c(";
                for (int j = 0; j < options.nCodes; j++) {
                    for (int k = 0; k < options.nCodes; k++) {
                        log << (j == 0 && k == 0 ? "" : ",") << strformat(" %.8g", transmat.eigeninv[j][k]);
                    }
                }
                log << ")" << std::endl;

                log << "# Transform into matrices and compute un-rotated vectors for profiles A and B" << std::endl;
                log << "codeFreq = matrix(codeFreq,nrow=20);" << std::endl;
                log << "eigeninv = matrix(eigeninv,nrow=20);" << std::endl;
                log << "unrotA = stat * (eigeninv %*% fA)" << std::endl;
                log << "unrotB = stat * (eigeninv %*% fB)" << std::endl;
                log << "# End of R block" << std::endl;

            }
            assert(lkAB > 0);
            if (site_likelihoods != nullptr)
                site_likelihoods[i] *= lkAB;
            lk *= lkAB;
            while (lk < Constants::LkUnderflow) {
                lk *= Constants::LkUnderflowInv;
                loglk -= Constants::LogLkUnderflow;
            }
            while (lk > Constants::LkUnderflowInv) {
                lk *= Constants::LkUnderflow;
                loglk += Constants::LogLkUnderflow;
            }
        }
    } else {
        assert(0);            /* illegal nCodes */
    }

    loglk += std::log(lk);
    options.debug.nLkCompute++;
    return (loglk);
}

AbsNeighbourJoining(double)::pairNegLogLk(double x, QuartetOpt &qo) {
    assert(qo.pair1 != nullptr && qo.pair2 != nullptr);
    qo.nEval++;
    double loglk = pairLogLk(*qo.pair1, *qo.pair2, x, /*site_lk*/nullptr);
    assert(loglk < 1e100);
    if (options.verbose > 5) {
        log << strformat("PairLogLk(%.4f) =  %.4f", x, loglk) << std::endl;
    }
    return -loglk;
}

AbsNeighbourJoining(void)::correctedPairDistances(Profile *_profiles[], int nProfiles, double distances[6]) {
    assert(_profiles != nullptr);
    assert(nProfiles > 1 && nProfiles <= 4);
    Besthit hit[6];

    for (int iHit = 0, i = 0; i < nProfiles; i++) {
        for (int j = i + 1; j < nProfiles; j++, iHit++) {
            profileDist(*_profiles[i], *_profiles[j], hit[iHit]);
            distances[iHit] = hit[iHit].dist;
        }
    }
    if (options.pseudoWeight > 0) {
        /* Estimate the prior distance */
        double dTop = 0;
        double dBottom = 0;
        for (int iHit = 0; iHit < (nProfiles * (nProfiles - 1)) / 2; iHit++) {
            dTop += hit[iHit].dist * hit[iHit].weight;
            dBottom += hit[iHit].weight;
        }
        double prior = (dBottom > 0.01) ? dTop / dBottom : 3.0;
        for (int iHit = 0; iHit < (nProfiles * (nProfiles - 1)) / 2; iHit++)
            distances[iHit] = (distances[iHit] * hit[iHit].weight + prior * options.pseudoWeight)
                              / (hit[iHit].weight + options.pseudoWeight);
    }
    if (options.logdist) {
        for (int iHit = 0; iHit < (nProfiles * (nProfiles - 1)) / 2; iHit++)
            distances[iHit] = logCorrect(distances[iHit]);
    }
}

AbsNeighbourJoining(void)::quartetConstraintPenalties(Profile *profiles4[4], double penalty[3]) {
    for (int i = 0; i < 3; i++) {
        penalty[i] = 0.0;
    }
    if (constraintSeqs.empty()) {
        return;
    }

    for (int64_t iC = 0; iC < (int64_t) constraintSeqs.size(); iC++) {
        double part[3];
        if (quartetConstraintPenaltiesPiece(profiles4, iC, /*OUT*/part)) {
            for (int i = 0; i < 3; i++) {
                penalty[i] += part[i];
            }

            if (options.verbose > 2
                && (std::fabs(part[ABvsCD] - part[ACvsBD]) > 0.001 || std::fabs(part[ABvsCD] - part[ADvsBC]) > 0.001)) {
                log << strformat(
                        "Constraint Penalties at %ld: ABvsCD %.3f ACvsBD %.3f ADvsBC %.3f %ld/%ld %ld/%ld %ld/%ld %ld/%ld",
                        iC, part[ABvsCD], part[ACvsBD], part[ADvsBC],
                        profiles4[0]->nOn[iC], profiles4[0]->nOff[iC],
                        profiles4[1]->nOn[iC], profiles4[1]->nOff[iC],
                        profiles4[2]->nOn[iC], profiles4[2]->nOff[iC],
                        profiles4[3]->nOn[iC], profiles4[3]->nOff[iC]) << std::endl;
            }
        }
    }
    if (options.verbose > 2) {
        log << strformat("Total Constraint Penalties: ABvsCD %.3f ACvsBD %.3f ADvsBC %.3f",
                         penalty[ABvsCD], penalty[ACvsBD], penalty[ADvsBC]) << std::endl;
    }
}

AbsNeighbourJoining(double)::pairConstraintDistance(int64_t nOn1, int64_t nOff1, int64_t nOn2, int64_t nOff2) {
    double f1 = nOn1 / (double) (nOn1 + nOff1);
    double f2 = nOn2 / (double) (nOn2 + nOff2);
    /* 1 - f1 * f2 - (1-f1)*(1-f2) = 1 - f1 * f2 - 1 + f1 + f2 - f1 * f2 */
    return (f1 + f2 - 2.0 * f1 * f2);
}

AbsNeighbourJoining(bool)::splitViolatesConstraint(Profile *profiles[4], int64_t iConstraint) {
    int codes[4]; /* 0 for off, 1 for on, -1 for split (quit if not constrained at all) */
    for (int64_t i = 0; i < 4; i++) {
        if (profiles[i]->nOn[iConstraint] + profiles[i]->nOff[iConstraint] == 0) {
            return false;
        } else if (profiles[i]->nOn[iConstraint] > 0 && profiles[i]->nOff[iConstraint] == 0) {
            codes[i] = 1;
        } else if (profiles[i]->nOn[iConstraint] == 0 && profiles[i]->nOff[iConstraint] > 0) {
            codes[i] = 0;
        } else {
            codes[i] = -1;
        }
    }
    int n0 = 0;
    int n1 = 0;
    for (int i = 0; i < 4; i++) {
        if (codes[i] == 0) {
            n0++;
        } else if (codes[i] == 1) {
            n1++;
        }
    }
    /* 3 on one side means no violation, even if other is code -1
       otherwise must have code != -1 and agreement on the split
     */
    if (n0 >= 3 || n1 >= 3) {
        return (false);
    }
    if (n0 == 2 && n1 == 2 && codes[0] == codes[1] && codes[2] == codes[3]) {
        return (false);
    }
    return true;
}

AbsNeighbourJoining(bool)::quartetConstraintPenaltiesPiece(Profile *profiles4[4], int64_t iC, double piece[3]) {
    int64_t nOn[4];
    int64_t nOff[4];

    int64_t nSplit = 0;
    int64_t nPlus = 0;
    int64_t nMinus = 0;

    for (int64_t i = 0; i < 4; i++) {
        nOn[i] = profiles4[i]->nOn[iC];
        nOff[i] = profiles4[i]->nOff[iC];
        if (nOn[i] + nOff[i] == 0)
            return (false);        /* ignore */
        else if (nOn[i] > 0 && nOff[i] > 0)
            nSplit++;
        else if (nOn[i] > 0)
            nPlus++;
        else
            nMinus++;
    }
    /* If just one of them is split or on the other side and the others all agree, also ignore */
    if (nPlus >= 3 || nMinus >= 3) {
        return false;
    }
    piece[ABvsCD] = options.constraintWeight
                    * (pairConstraintDistance(nOn[0], nOff[0], nOn[1], nOff[1])
                       + pairConstraintDistance(nOn[2], nOff[2], nOn[3], nOff[3]));
    piece[ACvsBD] = options.constraintWeight
                    * (pairConstraintDistance(nOn[0], nOff[0], nOn[2], nOff[2])
                       + pairConstraintDistance(nOn[1], nOff[1], nOn[3], nOff[3]));
    piece[ADvsBC] = options.constraintWeight
                    * (pairConstraintDistance(nOn[0], nOff[0], nOn[3], nOff[3])
                       + pairConstraintDistance(nOn[2], nOff[2], nOn[1], nOff[1]));
    return true;
}

AbsNeighbourJoining(void)::seqDist(char codes1[], char codes2[], Besthit &hit) {
    double top = 0;        /* summed over positions */
    int64_t nUse = 0;
    if (!distanceMatrix) {
        int nDiff = 0;
        for (int64_t i = 0; i < nPos; i++) {
            if (codes1[i] != NOCODE && codes2[i] != NOCODE) {
                nUse++;
                if (codes1[i] != codes2[i]) nDiff++;
            }
        }
        top = (double) nDiff;
    } else {
        for (int64_t i = 0; i < nPos; i++) {
            if (codes1[i] != NOCODE && codes2[i] != NOCODE) {
                nUse++;
                top += distanceMatrix.distances[(int64_t) codes1[i]][(int64_t) codes2[i]];
            }
        }
    }
    hit.weight = (double) nUse;
    hit.dist = nUse > 0 ? top / (double) nUse : 1.0;
    options.debug.seqOps++;
}

AbsNeighbourJoining(bool)::updateBestHit(int64_t nActive, Besthit &hit, bool bUpdateDist) {
    int64_t i = activeAncestor(hit.i);
    int64_t j = activeAncestor(hit.j);
    if (i < 0 || j < 0 || i == j) {
        hit.i = -1;
        hit.j = -1;
        hit.weight = 0;
        hit.dist = (numeric_t) 1e20;
        hit.criterion = (numeric_t) 1e20;
        return false;
    }
    if (i != hit.i || j != hit.j) {
        hit.i = i;
        hit.j = j;
        if (bUpdateDist) {
            setDistCriterion(nActive, hit);
        } else {
            hit.dist = (numeric_t) - 1e20;
            hit.criterion = (numeric_t) 1e20;
        }
    }
    return true;
}

AbsNeighbourJoining(double)::MLQuartetOptimize(Profile &pA, Profile &pB, Profile &pC, Profile &pD,
                                               double branch_lengths[5], bool *pStarTest, double *site_likelihoods) {
    double start_length[5];
    for (int j = 0; j < 5; j++) {
        start_length[j] = branch_lengths[j];
        if (branch_lengths[j] < options.MLMinBranchLength) {
            branch_lengths[j] = options.MLMinBranchLength;
        }
    }

    QuartetOpt qopt = {/*nEval*/0, /*pair1*/nullptr, /*pair2*/nullptr};
    double f2x, negloglk;

    if (pStarTest != nullptr) {
        *pStarTest = false;
    }

    /* First optimize internal branch, then branch to A, B, C, D, in turn
       May use star test to quit after internal branch
     */
    Profile AB(nPos, /*nConstraints*/0);
    Profile CD(nPos, /*nConstraints*/0);

    posteriorProfile(AB, pA, pB, branch_lengths[LEN_A], branch_lengths[LEN_B]);
    posteriorProfile(CD, pC, pD, branch_lengths[LEN_C], branch_lengths[LEN_D]);
    qopt.pair1 = &AB;
    qopt.pair2 = &CD;
    branch_lengths[LEN_I] = onedimenmin(
            /*xmin*/options.MLMinBranchLength,
            /*xguess*/branch_lengths[LEN_I],
            /*xmax*/6.0,
            /*func*/[this](double x, QuartetOpt &qo) { return pairNegLogLk(x, qo); },
            /*data*/qopt,
            /*ftol*/options.MLFTolBranchLength,
            /*atol*/options.MLMinBranchLengthTolerance,
            /*OUT*/negloglk,
            /*OUT*/f2x);

    if (pStarTest != nullptr) {
        assert(site_likelihoods == nullptr);
        double loglkStar = -pairNegLogLk(options.MLMinBranchLength, qopt);
        if (loglkStar < -negloglk - Constants::closeLogLkLimit) {
            *pStarTest = true;
            double off = pairLogLk(pA, pB, branch_lengths[LEN_A] + branch_lengths[LEN_B],/*site_lk*/nullptr)
                         + pairLogLk(pC, pD, branch_lengths[LEN_C] + branch_lengths[LEN_D], /*site_lk*/nullptr);
            return -negloglk + off;
        }
    }
    {
        Profile BCD(nPos, /*nConstraints*/0);
        posteriorProfile(BCD, pB, CD, branch_lengths[LEN_B], branch_lengths[LEN_I]);
        qopt.pair1 = &pA;
        qopt.pair2 = &BCD;
        branch_lengths[LEN_A] = onedimenmin(
                /*xmin*/options.MLMinBranchLength,
                /*xguess*/branch_lengths[LEN_A],
                /*xmax*/6.0,
                /*func*/[this](double x, QuartetOpt &qo) { return pairNegLogLk(x, qo); },
                /*data*/qopt,
                /*ftol*/options.MLFTolBranchLength,
                /*atol*/options.MLMinBranchLengthTolerance,
                /*OUT*/negloglk,
                /*OUT*/f2x);
    }
    {
        Profile ACD(nPos, /*nConstraints*/0);
        posteriorProfile(ACD, pA, CD, branch_lengths[LEN_A], branch_lengths[LEN_I]);
        qopt.pair1 = &pB;
        qopt.pair2 = &ACD;
        branch_lengths[LEN_B] = onedimenmin(
                /*xmin*/options.MLMinBranchLength,
                /*xguess*/branch_lengths[LEN_B],
                /*xmax*/6.0,
                /*func*/[this](double x, QuartetOpt &qo) { return pairNegLogLk(x, qo); },
                /*data*/qopt,
                /*ftol*/options.MLFTolBranchLength,
                /*atol*/options.MLMinBranchLengthTolerance,
                /*OUT*/negloglk,
                /*OUT*/f2x);
    }
    posteriorProfile(AB, pA, pB, branch_lengths[LEN_A], branch_lengths[LEN_B]);
    {
        Profile ABD(nPos, /*nConstraints*/0);
        posteriorProfile(ABD, AB, pD, branch_lengths[LEN_I], branch_lengths[LEN_D]);
        qopt.pair1 = &pC;
        qopt.pair2 = &ABD;
        branch_lengths[LEN_C] = onedimenmin(
                /*xmin*/options.MLMinBranchLength,
                /*xguess*/branch_lengths[LEN_C],
                /*xmax*/6.0,
                /*func*/[this](double x, QuartetOpt &qo) { return pairNegLogLk(x, qo); },
                /*data*/qopt,
                /*ftol*/options.MLFTolBranchLength,
                /*atol*/options.MLMinBranchLengthTolerance,
                /*OUT*/negloglk,
                /*OUT*/f2x);
    }

    Profile ABC(nPos, /*nConstraints*/0);
    posteriorProfile(ABC, AB, pC, branch_lengths[LEN_I], branch_lengths[LEN_C]);
    qopt.pair1 = &pD;
    qopt.pair2 = &ABC;
    branch_lengths[LEN_D] = onedimenmin(
            /*xmin*/options.MLMinBranchLength,
            /*xguess*/branch_lengths[LEN_D],
            /*xmax*/6.0,
            /*func*/[this](double x, QuartetOpt &qo) { return pairNegLogLk(x, qo); },
            /*data*/qopt,
            /*ftol*/options.MLFTolBranchLength,
            /*atol*/options.MLMinBranchLengthTolerance,
            /*OUT*/negloglk,
            /*OUT*/f2x);

    /* Compute the total quartet likelihood
       PairLogLk(ABC,D) + PairLogLk(AB,C) + PairLogLk(A,B)
     */
    double loglkABCvsD = -negloglk;
    if (site_likelihoods) {
        for (int64_t j = 0; j < nPos; j++) {
            site_likelihoods[j] = 1.0;
        }
        pairLogLk(ABC, pD, branch_lengths[LEN_D],/*IN/OUT*/site_likelihoods);
    }
    double quartetloglk = loglkABCvsD
                          + pairLogLk(AB, pC, branch_lengths[LEN_I] + branch_lengths[LEN_C],/*IN/OUT*/site_likelihoods)
                          + pairLogLk(pA, pB, branch_lengths[LEN_A] + branch_lengths[LEN_B],/*IN/OUT*/site_likelihoods);

    if (options.verbose > 3) {
        double loglkStart = MLQuartetLogLk(pA, pB, pC, pD, start_length, /*site_lk*/nullptr);
        log << strformat("Optimize loglk from %.5f to %.5f eval %ld lengths from\n"
                         "   %.5f %.5f %.5f %.5f %.5f to\n"
                         "   %.5f %.5f %.5f %.5f %.5f",
                         loglkStart, quartetloglk, qopt.nEval,
                         start_length[0], start_length[1], start_length[2], start_length[3], start_length[4],
                         branch_lengths[0], branch_lengths[1], branch_lengths[2], branch_lengths[3], branch_lengths[4])
            << std::endl;
    }
    return quartetloglk;
}

AbsNeighbourJoining(double)::MLPairOptimize(Profile &pA, Profile &pB, double *branch_length) {
    QuartetOpt qopt = {/*nEval*/0, /*pair1*/&pA, /*pair2*/&pB};
    double f2x, negloglk;
    *branch_length = onedimenmin(/*xmin*/options.MLMinBranchLength,
            /*xguess*/*branch_length,
            /*xmax*/6.0,
            /*func*/[this](double x, QuartetOpt &qo) { return pairNegLogLk(x, qo); },
            /*data*/qopt,
            /*ftol*/options.MLFTolBranchLength,
            /*atol*/options.MLMinBranchLengthTolerance,
            /*OUT*/negloglk,
            /*OUT*/f2x);
    return -negloglk;        /* the log likelihood */
}

AbsNeighbourJoining(int64_t)::findSPRSteps(int64_t nodeMove, int64_t nodeAround, std::unique_ptr<Profile> upProfiles[],
                                           SprStep steps[], bool bFirstAC) {
    int64_t iStep;
    for (iStep = 0; iStep < options.maxSPRLength; iStep++) {
        if (child[nodeAround].nChild != 2) {
            break;            /* no further to go */
        }

        /* Consider the NNIs around nodeAround */
        Profile *profiles[4];
        int64_t nodeABCD[4];
        setupABCD(nodeAround, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, /*useML*/false);
        double criteria[3];
        chooseNNI(profiles, criteria);

        /* Do & save the swap */
        auto &step = steps[iStep];
        if (iStep == 0 ? bFirstAC : criteria[ACvsBD] < criteria[ADvsBC]) {
            /* swap B & C to put AC together */
            step.deltaLength = criteria[ACvsBD] - criteria[ABvsCD];
            step.nodes[0] = nodeABCD[1];
            step.nodes[1] = nodeABCD[2];
        } else {
            /* swap AC to put AD together */
            step.deltaLength = criteria[ADvsBC] - criteria[ABvsCD];
            step.nodes[0] = nodeABCD[0];
            step.nodes[1] = nodeABCD[2];
        }

        if (options.verbose > 3) {
            log << strformat("SPR chain step %ld for %ld around %ld swap %ld %ld deltaLen %.5f",
                             iStep + 1, nodeAround, nodeMove, step.nodes[0], step.nodes[1], step.deltaLength)
                << std::endl;
            if (options.verbose > 4)
                printNJInternal(log, /*useLen*/false);
        }
        replaceChild(nodeAround, step.nodes[0], step.nodes[1]);
        replaceChild(parent[nodeAround], step.nodes[1], step.nodes[0]);
        updateForNNI(nodeAround, /*IN/OUT*/upProfiles, /*useML*/false);

        /* set the new nodeAround -- either parent(nodeMove) or sibling(nodeMove) --
           so that it different from current nodeAround
         */
        int64_t newAround[2] = {parent[nodeMove], sibling(nodeMove)};
        if (parent[nodeMove] == root) {
            rootSiblings(nodeMove, /*OUT*/newAround);
        }
        assert(newAround[0] == nodeAround || newAround[1] == nodeAround);
        assert(newAround[0] != newAround[1]);
        nodeAround = newAround[newAround[0] == nodeAround ? 1 : 0];
    }
    return iStep;
}


AbsNeighbourJoining(void)::unwindSPRStep(SprStep &step, std::unique_ptr<Profile> upProfiles[]) {
    int64_t parents[2];
    for (int i = 0; i < 2; i++) {
        assert(step.nodes[i] >= 0 && step.nodes[i] < maxnodes);
        parents[i] = parent[step.nodes[i]];
        assert(parents[i] >= 0);
    }
    assert(parents[0] != parents[1]);
    replaceChild(parents[0], step.nodes[0], step.nodes[1]);
    replaceChild(parents[1], step.nodes[1], step.nodes[0]);
    int iYounger = 0;
    if (parent[parents[0]] == parents[1]) {
        iYounger = 0;
    } else {
        assert(parent[parents[1]] == parents[0]);
        iYounger = 1;
    }
    updateForNNI(parents[iYounger], /*IN/OUT*/upProfiles, /*useML*/false);
}

/* Update the profile of node and its ancestor, and delete nearby out-profiles */
AbsNeighbourJoining(void)::updateForNNI(int64_t node, std::unique_ptr<Profile> upProfiles[], bool useML) {
    if (options.slow) {
        /* exhaustive update */
        for (int64_t i = 0; i < maxnodes; i++) {
            upProfiles[i].reset();
        }

        /* update profiles back to root */
        for (int64_t ancestor = node; ancestor >= 0; ancestor = parent[ancestor]) {
            /* parallel SPR cannot modify outside its range. */
            if (!lockedNodes.empty() && lockedNodes[ancestor]) {
                break;
            }
            recomputeProfile(upProfiles, ancestor, useML);
        }

        /* remove any up-profiles made while doing that*/
        for (int64_t i = 0; i < maxnodes; i++) {
            upProfiles[i].reset();
        }
    } else {
        /* if fast, only update around self
           note that upProfile(parent) is still OK after an NNI, but
           up-profiles of uncles may not be
        */
        upProfiles[node].reset();
        for (int64_t i = 0; i < child[node].nChild; i++) {
            upProfiles[child[node].child[i]].reset();
        }
        assert(node != root);
        int64_t iparent = parent[node];
        int64_t neighbors[2] = {iparent, sibling(node)};
        if (iparent == root) {
            rootSiblings(node, /*OUT*/neighbors);
        }
        upProfiles[neighbors[0]].reset();
        upProfiles[neighbors[1]].reset();

        int64_t uncle = sibling(iparent);
        if (uncle >= 0) {
            upProfiles[uncle].reset();
        }
        recomputeProfile(upProfiles, node, useML);
        recomputeProfile(upProfiles, iparent, useML);
    }
}

/* Resets the children entry of parent and also the parent entry of newchild */
AbsNeighbourJoining(void)::replaceChild(int64_t _parent, int64_t oldchild, int64_t newchild) {
    parent[newchild] = _parent;

    for (int64_t iChild = 0; iChild < child[_parent].nChild; iChild++) {
        if (child[_parent].child[iChild] == oldchild) {
            child[_parent].child[iChild] = newchild;
            return;
        }
    }
    assert(0);
}

AbsNeighbourJoining(void)::setupABCD(int64_t node, Profile *profiles4[4], std::unique_ptr<Profile> upProfiles[],
                                     int64_t nodeABCD[4], bool useML) {
    int64_t iparent = parent[node];
    assert(iparent >= 0);
    assert(child[node].nChild == 2);
    nodeABCD[0] = child[node].child[0]; /*A*/
    nodeABCD[1] = child[node].child[1]; /*B*/

    Profile *profile = nullptr;
    if (iparent == root) {
        int64_t sibs[2];
        rootSiblings(node, /*OUT*/sibs);
        nodeABCD[2] = sibs[0];
        nodeABCD[3] = sibs[1];
        if (profiles4 == nullptr) {
            return;
        }
        profile = &profiles[sibs[1]];
    } else {
        nodeABCD[2] = sibling(node);
        assert(nodeABCD[2] >= 0);
        nodeABCD[3] = iparent;
        if (profiles4 == nullptr) {
            return;
        }
        profile = getUpProfile(upProfiles, iparent, useML);
    }
    assert(upProfiles);
    for (int i = 0; i < 3; i++) {
        profiles4[i] = &profiles[nodeABCD[i]];
    }
    profiles4[3] = profile;
}

AbsNeighbourJoining(int64_t)::sibling(int64_t node) {
    int64_t iparent = parent[node];
    if (iparent < 0 || iparent == root) {
        return -1;
    }

    for (int64_t iChild = 0; iChild < child[iparent].nChild; iChild++) {
        if (child[iparent].child[iChild] != node) {
            return (child[iparent].child[iChild]);
        }
    }
    assert(0);
    return -1;
}

AbsNeighbourJoining(void)::rootSiblings(int64_t node, /*OUT*/int64_t sibs[2]) {
    assert(parent[node] == root);
    assert(child[root].nChild == 3);

    int64_t nSibs = 0;
    for (int64_t iChild = 0; iChild < child[root].nChild; iChild++) {
        int64_t ichild = child[root].child[iChild];
        if (ichild != node) {
            assert(nSibs <= 2);
            sibs[nSibs++] = ichild;
        }
    }
}

AbsNeighbourJoining(void)::pSameVector(double length, std::vector<double> &pSame) {
    pSame.resize(rates.rates.size()); /* Value is 1 or 0 so vector operations can not be used*/
    for (int64_t iRate = 0; iRate < (int64_t) rates.rates.size(); iRate++) {
        pSame[iRate] = 0.25 + 0.75 * std::exp((-4.0 / 3.0) * std::abs(length * rates.rates[iRate]));
    }
}


AbsNeighbourJoining(void)::pDiffVector(std::vector<double> &pSame, std::vector<double> &pDiff) {
    pDiff.resize(rates.rates.size());
    for (int64_t iRate = 0; iRate < (int64_t) rates.rates.size(); iRate++) {
        pDiff[iRate] = (1.0 - pSame[iRate]) / 3.0;
    }
}

AbsNeighbourJoining(void)::
expEigenRates(double length, std::vector<numeric_t, typename op_t::Allocator> &expeigenRates) {
    expeigenRates.resize(nCodeSize * rates.rates.size());
    for (int64_t iRate = 0; iRate < (int64_t) rates.rates.size(); iRate++) {
        double relLen = length * rates.rates[iRate];
        /* very short branch lengths lead to numerical problems so prevent them */
        if (relLen < options.MLMinRelBranchLength) {
            relLen = options.MLMinRelBranchLength;
        }
#ifndef NDEBUG
        for (int64_t j = 0; j < options.nCodes; j++) {
            expeigenRates[iRate * nCodeSize + j] = std::exp((double) relLen * transmat.eigenval[j]);
        }
#else
        operations.vector_multiply_by(transmat.eigenval, relLen, options.nCodes, &expeigenRates[iRate * nCodeSize]);
        operations.fastexp(&expeigenRates[iRate * nCodeSize], options.nCodes, options.fastexp);
#endif
    }
}

AbsNeighbourJoining(inline Precision*)::getFreq(Profile &p, int64_t &i, int64_t &ivector) {
    return p.weights[i] > 0 && p.codes[i] == NOCODE ? &p.vectors[nCodeSize * (ivector++)] : nullptr;
}

AbsNeighbourJoining(int64_t)::nGaps(int64_t i) {
    assert(i < (int64_t) nSeqs);
    int nGaps = 0;
    for (int64_t p = 0; p < nPos; p++) {
        if (profiles[i].codes[p] == NOCODE) {
            nGaps++;
        }
    }
    return nGaps;
}

AbsNeighbourJoining(void)::setProfile(int64_t node, double weight1) {
    Children &c = child[node];
    assert(c.nChild == 2);
    assert((int64_t) profiles.size() > c.child[0]);
    assert((int64_t) profiles.size() > c.child[1]);
    averageProfile(profiles[node], profiles[c.child[0]], profiles[c.child[1]], weight1);
}

AbsNeighbourJoining(void)::averageProfile(Profile &out, Profile &profile1, Profile &profile2, double weight1) {
    averageProfile(out, profile1, profile2, weight1, distanceMatrix);
}

AbsNeighbourJoining(void)::averageProfile(Profile &out, Profile &profile1, Profile &profile2, double bionjWeight,
                                          DistanceMatrix <Precision, op_t::ALIGNMENT> &dmat) {
    if (bionjWeight < 0) {
        bionjWeight = 0.5;
    }
    out.reset();

    for (int64_t i = 0; i < nPos; i++) {
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
                out.nVectors++;
            }
        }
    }

    /* Allocate and set the vectors */
    out.setVectorSize(out.nVectors * nCodeSize, 0);

    options.debug.nProfileFreqAlloc += out.nVectors;
    options.debug.nProfileFreqAvoid += nPos - out.nVectors;
    int64_t iFreqOut = 0;
    int64_t iFreq1 = 0;
    int64_t iFreq2 = 0;
    for (int64_t i = 0; i < nPos; i++) {
        numeric_t *f = getFreq(out, i, iFreqOut);
        numeric_t *f1 = getFreq(profile1, i, iFreq1);
        numeric_t *f2 = getFreq(profile2, i, iFreq2);
        if (f != nullptr) {
            if (profile1.weights[i] > 0) {
                addToFreq(f, profile1.weights[i] * bionjWeight, profile1.codes[i], f1, dmat);
            }
            if (profile2.weights[i] > 0) {
                addToFreq(f, profile2.weights[i] * (1.0 - bionjWeight), profile2.codes[i], f2, dmat);
            }
            normalizeFreq(f, dmat);
        } /* end if computing f */
        if (options.verbose > 10 && i < 5) {
            log << strformat("Average profiles: pos %ld in-w1 %f in-w2 %f bionjWeight %f to weight %f code %d",
                             i, profile1.weights[i], profile2.weights[i], bionjWeight, out.weights[i], out.codes[i])
                << std::endl;
            if (f != nullptr) {
                for (int k = 0; k < options.nCodes; k++) {
                    log << strformat("\t%c:%f", options.codesString[k], f ? f[k] : -1.0);
                }
                log << std::endl;
            }
        }
    } /* end loop over positions */
    assert(iFreq1 == profile1.nVectors);
    assert(iFreq2 == profile2.nVectors);
    assert(iFreqOut == out.nVectors);

    /* compute total constraints */
    for (int64_t i = 0; i < (int64_t) constraintSeqs.size(); i++) {
        out.nOn[i] = profile1.nOn[i] + profile2.nOn[i];
        out.nOff[i] = profile1.nOff[i] + profile2.nOff[i];
    }
    options.debug.profileAvgOps++;
}

AbsNeighbourJoining(void)::
posteriorProfile(Profile &out, Profile &p1, Profile &p2, double len1, double len2) {
    if (len1 < options.MLMinBranchLength) {
        len1 = options.MLMinBranchLength;
    }
    if (len2 < options.MLMinBranchLength) {
        len2 = options.MLMinBranchLength;
    }

    for (int64_t i = 0; i < nPos; i++) {
        out.codes[i] = NOCODE;
        out.weights[i] = 1.0;
    }

    out.nVectors = nPos;
    out.setVectorSize(out.nVectors * nCodeSize, 0);

    int64_t iFreqOut = 0;
    int64_t iFreq1 = 0;
    int64_t iFreq2 = 0;
    std::vector<numeric_t, typename op_t::Allocator> expeigenRates1, expeigenRates2;

    if (transmat) {
        expEigenRates(len1, expeigenRates1);
        expEigenRates(len2, expeigenRates2);
    }

    if (!transmat) {    /* Jukes-Cantor */
        assert(options.nCodes == 4);

        std::vector<double> PSame1, PDiff1, PSame2, PDiff2;

        pSameVector(len1, PSame1);
        pDiffVector(PSame1, PDiff1);
        pSameVector(len2, PSame2);
        pDiffVector(PSame2, PDiff2);

        numeric_t mix1[4], mix2[4];

        for (int64_t i = 0; i < nPos; i++) {
            int64_t iRate = rates.ratecat[i];
            double w1 = p1.weights[i];
            double w2 = p2.weights[i];
            int code1 = p1.codes[i];
            int code2 = p2.codes[i];
            numeric_t *f1 = getFreq(p1, i,/*IN/OUT*/iFreq1);
            numeric_t *f2 = getFreq(p2, i,/*IN/OUT*/iFreq2);

            /* First try to store a simple profile */
            if (f1 == nullptr && f2 == nullptr) {
                if (code1 == NOCODE && code2 == NOCODE) {
                    out.codes[i] = NOCODE;
                    out.weights[i] = 0.0;
                    continue;
                } else if (code1 == NOCODE) {
                    /* Posterior(parent | character & gap, len1, len2) = Posterior(parent | character, len1)
                       = PSame() for matching characters and 1-PSame() for the rest
                       = (pSame - pDiff) * character + (1-(pSame-pDiff)) * gap
                    */
                    out.codes[i] = (char) code2;
                    out.weights[i] = w2 * (PSame2[iRate] - PDiff2[iRate]);
                    continue;
                } else if (code2 == NOCODE) {
                    out.codes[i] = (char) code1;
                    out.weights[i] = w1 * (PSame1[iRate] - PDiff1[iRate]);
                    continue;
                } else if (code1 == (char) code2) {
                    out.codes[i] = (char) code1;
                    double f12code = (w1 * PSame1[iRate] + (1 - w1) * 0.25) * (w2 * PSame2[iRate] + (1 - w2) * 0.25);
                    double f12other = (w1 * PDiff1[iRate] + (1 - w1) * 0.25) * (w2 * PDiff2[iRate] + (1 - w2) * 0.25);
                    /* posterior probability of code1/code2 after scaling */
                    double pcode = f12code / (f12code + 3 * f12other);
                    /* Now f = w * (code ? 1 : 0) + (1-w) * 0.25, so to get pcode we need
                       fcode = 1/4 + w1*3/4 or w = (f-1/4)*4/3
                     */
                    out.weights[i] = (pcode - 0.25) * 4.0 / 3.0;
                    /* This can be zero because of numerical problems, I think */
                    if (out.weights[i] < 1e-6) {
                        if (options.verbose > 1) {
                            log << strformat(
                                    "Replaced weight %f with %f from w1 %f w2 %f PSame %f %f f12code %f f12other %f",
                                    out.weights[i], 1e-6,
                                    w1, w2,
                                    PSame1[iRate], PSame2[iRate],
                                    f12code, f12other) << std::endl;
                        }
                        out.weights[i] = (numeric_t) 1e-6;
                    }
                    continue;
                }
            }
            /* if we did not compute a simple profile, then do the full computation and
               store the full vector
            */
            if (f1 == nullptr) {
                for (int j = 0; j < 4; j++) {
                    mix1[j] = (1 - w1) * 0.25;
                }
                if (code1 != NOCODE)
                    mix1[code1] += w1;
                f1 = mix1;
            }
            if (f2 == nullptr) {
                for (int j = 0; j < 4; j++) {
                    mix2[j] = (1 - w2) * 0.25;
                }
                if (code2 != NOCODE) {
                    mix2[code2] += w2;
                }
                f2 = mix2;
            }
            out.codes[i] = NOCODE;
            out.weights[i] = 1.0;
            numeric_t *f = getFreq(out, i,/*IN/OUT*/iFreqOut);
            double lkAB = 0;
            for (int j = 0; j < 4; j++) {
                f[j] = (f1[j] * PSame1[iRate] + (1.0 - f1[j]) * PDiff1[iRate])
                       * (f2[j] * PSame2[iRate] + (1.0 - f2[j]) * PDiff2[iRate]);
                lkAB += f[j];
            }
            double lkABInv = 1.0 / lkAB;
            for (int j = 0; j < 4; j++) {
                f[j] *= lkABInv;
            }
        }
    } else if (options.nCodes == 4) {    /* matrix model on nucleotides */
        numeric_t *fGap = &transmat.codeFreq[NOCODE][0];
        numeric_t f1mix[4], f2mix[4];

        for (int64_t i = 0; i < nPos; i++) {
            if (p1.codes[i] == NOCODE && p2.codes[i] == NOCODE && p1.weights[i] == 0 && p2.weights[i] == 0) {
                /* aligning gap with gap -- just output a gap
                   out->codes[i] is already set to NOCODE so need not set that */
                out.weights[i] = 0;
                continue;
            }
            int64_t iRate = rates.ratecat[i];
            numeric_t *expeigen1 = &expeigenRates1[iRate * nCodeSize];
            numeric_t *expeigen2 = &expeigenRates2[iRate * nCodeSize];
            numeric_t *f1 = getFreq(p1, i,/*IN/OUT*/iFreq1);
            numeric_t *f2 = getFreq(p2, i,/*IN/OUT*/iFreq2);
            numeric_t *fOut = getFreq(out, i,/*IN/OUT*/iFreqOut);
            assert(fOut != nullptr);

            if (f1 == nullptr) {
                f1 = &transmat.codeFreq[(int) p1.codes[i]][0]; /* codeFreq includes an entry for NOCODE */
                double w = p1.weights[i];
                if (w > 0.0 && w < 1.0) {
                    for (int j = 0; j < 4; j++) {
                        f1mix[j] = w * f1[j] + (1.0 - w) * fGap[j];
                    }
                    f1 = f1mix;
                }
            }
            if (f2 == nullptr) {
                f2 = &transmat.codeFreq[(int) p2.codes[i]][0];
                double w = p2.weights[i];
                if (w > 0.0 && w < 1.0) {
                    for (int j = 0; j < 4; j++) {
                        f2mix[j] = w * f2[j] + (1.0 - w) * fGap[j];
                    }
                    f2 = f2mix;
                }
            }

            numeric_t fMult1[4];  /* rotated1 * expeigen1 */
            numeric_t fMult2[4]; /* rotated2 * expeigen2 */

            for (int j = 0; j < 4; j++) {
                fMult1[j] = f1[j] * expeigen1[j];
                fMult2[j] = f2[j] * expeigen2[j];
            }

            alignas(op_t::ALIGNMENT) numeric_t fPost[4];  /* in  unrotated space */
            for (int j = 0; j < 4; j++) {
                double out1 = 0;
                double out2 = 0;
                for (int k = 0; k < 4; k++) {
                    out1 += fMult1[k] * transmat.codeFreq[j][k];
                    out2 += fMult2[k] * transmat.codeFreq[j][k];
                }
                fPost[j] = out1 * out2 * transmat.statinv[j];

            }
            double fPostTot = 0;
            for (int j = 0; j < 4; j++) {
                fPostTot += fPost[j];
            }
            assert(fPostTot > options.fPostTotalTolerance);
            double fPostInv = 1.0 / fPostTot;

            for (int j = 0; j < 4; j++) {
                fPost[j] *= fPostInv;
            }

            /* and finally, divide by stat again & rotate to give the new frequencies */
            operations.matrixt_by_vector4(transmat.eigeninvT, fPost, fOut);
        }  /* end loop over position i */
    } else if (options.nCodes == 20) {    /* matrix model on amino acids */
        numeric_t *fGap = &transmat.codeFreq[NOCODE][0];
        alignas(op_t::ALIGNMENT) numeric_t f1mix[20];
        alignas(op_t::ALIGNMENT) numeric_t f2mix[20];

        for (int64_t i = 0; i < nPos; i++) {
            if (p1.codes[i] == NOCODE && p2.codes[i] == NOCODE && p1.weights[i] == 0 && p2.weights[i] == 0) {
                /* aligning gap with gap -- just output a gap
                   out->codes[i] is already set to NOCODE so need not set that */
                out.weights[i] = 0;
                continue;
            }
            int64_t iRate = rates.ratecat[i];
            numeric_t *expeigen1 = &expeigenRates1[iRate * nCodeSize];
            numeric_t *expeigen2 = &expeigenRates2[iRate * nCodeSize];
            numeric_t *f1 = getFreq(p1, i, iFreq1);
            numeric_t *f2 = getFreq(p2, i, iFreq2);
            numeric_t *fOut = getFreq(out, i, iFreqOut);
            assert(fOut != nullptr);

            if (f1 == nullptr) {
                f1 = &transmat.codeFreq[(int) p1.codes[i]][0]; /* codeFreq includes an entry for NOCODE */
                double w = p1.weights[i];
                if (w > 0.0 && w < 1.0) {
                    for (int j = 0; j < 20; j++) {
                        f1mix[j] = w * f1[j] + (1.0 - w) * fGap[j];
                    }
                    f1 = f1mix;
                }
            }
            if (f2 == nullptr) {
                f2 = &transmat.codeFreq[(int) p2.codes[i]][0];
                double w = p2.weights[i];
                if (w > 0.0 && w < 1.0) {
                    for (int j = 0; j < 20; j++) {
                        f2mix[j] = w * f2[j] + (1.0 - w) * fGap[j];
                    }
                    f2 = f2mix;
                }
            }
            alignas(op_t::ALIGNMENT) numeric_t fMult1[20];    /* rotated1 * expeigen1 */
            alignas(op_t::ALIGNMENT) numeric_t fMult2[20];    /* rotated2 * expeigen2 */
            operations.vector_multiply(f1, expeigen1, 20, /*OUT*/fMult1);
            operations.vector_multiply(f2, expeigen2, 20, /*OUT*/fMult2);
            alignas(op_t::ALIGNMENT) numeric_t fPost[20];        /* in  unrotated space */
            for (int j = 0; j < 20; j++) {
                numeric_t value = operations.vector_dot_product_rot(fMult1, fMult2, &transmat.codeFreq[j][0], 20)
                                  * transmat.statinv[j];
                /* Added this logic try to avoid rare numerical problems */
                fPost[j] = value >= 0 ? value : 0;
            }
            double fPostTot = operations.vector_sum(fPost, 20);
            assert(fPostTot > options.fPostTotalTolerance);
            double fPostInv = 1.0 / fPostTot;
            operations.vector_multiply_by(fPost, fPostInv, 20, fPost);
            int64_t ch = -1;        /* the dominant character, if any */
            if (!options.exactML) {
                for (int j = 0; j < 20; j++) {
                    if (fPost[j] >= Constants::approxMLminf) {
                        ch = j;
                        break;
                    }
                }
            }

            /* now, see if we can use the approximation
           fPost ~= (1 or 0) * w + nearP * (1-w)
           to avoid rotating */
            double w = 0;
            if (ch >= 0) {
                w = (fPost[ch] - transmat.nearP[ch][ch]) / (1.0 - transmat.nearP[ch][ch]);
                for (int j = 0; j < 20; j++) {
                    if (j != ch) {
                        double fRough = (1.0 - w) * transmat.nearP[ch][j];
                        if (fRough < fPost[j] * Constants::approxMLminratio) {
                            ch = -1;        /* give up on the approximation */
                            break;
                        }
                    }
                }
            }
            if (ch >= 0) {
                options.debug.nAAPosteriorRough++;
                double wInvStat = w * transmat.statinv[ch];
                for (int j = 0; j < 20; j++) {
                    fOut[j] = wInvStat * transmat.codeFreq[ch][j] + (1.0 - w) * transmat.nearFreq[ch][j];
                }
            } else {
                /* and finally, divide by stat again & rotate to give the new frequencies */
                options.debug.nAAPosteriorExact++;
                for (int j = 0; j < 20; j++) {
                    fOut[j] = operations.vector_multiply_sum(fPost, &transmat.eigeninv[j][0], 20);
                }
            }
        } /* end loop over position i */
    } else {
        assert(0);            /* illegal nCodes */
    }

    /* Reallocate out->vectors to be the right size */
    out.nVectors = iFreqOut;
    out.resizeVector(iFreqOut * nCodeSize); /* try to save space */

    options.debug.nProfileFreqAlloc += out.nVectors;
    options.debug.nProfileFreqAvoid += nPos - out.nVectors;

    /* compute total constraints */
    for (int64_t i = 0; i < (int64_t) constraintSeqs.size(); i++) {
        out.nOn[i] = p1.nOn[i] + p2.nOn[i];
        out.nOff[i] = p1.nOff[i] + p2.nOff[i];
    }
    options.debug.nPosteriorCompute++;
}

AbsNeighbourJoining(void)::readTree(Uniquify &unique, HashTable &hashnames, std::istream &fpInTree) {
    assert(nSeqs == unique.uniqueSeq.size());
    /* First, do a preliminary parse of the tree to with non-unique leaves ignored
       We need to store this separately from NJ because it may have too many internal nodes
       (matching sequences show up once in the NJ but could be in multiple places in the tree)
       Will use iUnique as the index of nodes, as in the NJ structure
    */
    int64_t maxnode = unique.alnToUniq.size();
    int64_t maxnodes = maxnode * 2;
    int64_t root = maxnode++;
    std::vector<int64_t> parent(maxnodes, -1);
    std::vector<Children> children(maxnodes);

    /* The stack is the current path to the root, with the root at the first (top) position */
    int64_t stack_size = 1;
    std::vector<int64_t> stack(maxnodes);
    stack[0] = root;
    int64_t nDown = 0;
    int64_t nUp = 0;

    std::string token;
    token.reserve(5000l);

    if (!readTreeToken(fpInTree, token) || token[0] != '(') {
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
                    int64_t newnode = maxnode++;
                    assert(newnode < maxnodes);
                    readTreeAddChild(stack[stack_size - 1], newnode, parent, children);
                    if (options.verbose > 5) {
                        log << strformat("Added internal child %ld of %ld, stack size increase to %ld",
                                         newnode, stack[stack_size - 1], stack_size + 1) << std::endl;
                    }
                    stack[stack_size++] = newnode;
                    assert(stack_size < maxnodes);
                }
                readTreeMaybeAddLeaf(stack[stack_size - 1], token, hashnames, unique, parent, children);
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
                        log << strformat("Up to nUp=%ld stack size %ld at %ld",
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
            readTreeMaybeAddLeaf(stack[stack_size - 1], token, hashnames, unique, parent, children);
        }
    }

    /* Verify that all sequences were seen */
    for (int64_t i = 0; i < (int64_t) unique.uniqueSeq.size(); i++) {
        if (parent[i] < 0) {
            throw std::invalid_argument(
                    strformat("Alignment sequence %ld (unique %ld) absent from input tree\n"
                              "The starting tree (the argument to -intree) must include all sequences in the alignment!",
                              unique.uniqueFirst[i], i));
        }
    }

    /* Simplify the tree -- remove all internal nodes with < 2 children
       Keep trying until no nodes get removed
    */
    int64_t nRemoved;
    do {
        nRemoved = 0;
        /* Here stack is the list of nodes we haven't visited yet while doing
           a tree traversal */
        stack_size = 1;
        stack[0] = root;
        while (stack_size > 0) {
            int64_t node = stack[--stack_size];
            if (node >= (int64_t) unique.uniqueSeq.size()) { /* internal node */
                if (children[node].nChild <= 1) {
                    if (node != root) {
                        readTreeRemove(parent, children, node);
                        nRemoved++;
                    } else if (node == root && children[node].nChild == 1) {
                        int64_t newroot = children[node].child[0];
                        parent[newroot] = -1;
                        children[root].nChild = 0;
                        nRemoved++;
                        if (options.verbose > 5) {
                            log << strformat("Changed root from %ld to %ld", root, newroot) << std::endl;
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
                            log << strformat("Added %ld to stack", stack[stack_size - 1]) << std::endl;
                        }
                    }
                }
            }
        }
    } while (nRemoved > 0);

    /* Simplify the root node to 3 children if it has 2 */
    if (children[root].nChild == 2) {
        for (int64_t i = 0; i < 2; i++) {
            int64_t ichild = children[root].child[i];
            assert(ichild >= 0 && ichild < maxnodes);
            if (children[ichild].nChild == 2) {
                readTreeRemove(parent, children, ichild); /* replace root -> child -> A,B with root->A,B */
                break;
            }
        }
    }

    for (int64_t i = 0; i < maxnodes; i++)
        if (options.verbose > 5) {
            log << strformat("Simplfied node %ld has parent %ld nchild %d", i, parent[i], children[i].nChild)
                << std::endl;
        }

    /* Map the remaining internal nodes to NJ nodes */
    std::vector<int64_t> map(maxnodes);
    for (int64_t i = 0; i < (int64_t) unique.uniqueSeq.size(); i++) {
        map[i] = i;
    }
    for (int64_t i = unique.uniqueSeq.size(); i < maxnodes; i++) {
        map[i] = -1;
    }
    stack_size = 1;
    stack[0] = root;
    while (stack_size > 0) {
        int64_t node = stack[--stack_size];
        if (node >= (int64_t) unique.uniqueSeq.size()) { /* internal node */
            assert(node == root || children[node].nChild > 1);
            map[node] = this->maxnode++;
            for (int64_t i = 0; i < children[node].nChild; i++) {
                assert(stack_size < maxnodes);
                stack[stack_size++] = children[node].child[i];
            }
        }
    }
    for (int64_t i = 0; i < maxnodes; i++)
        if (options.verbose > 5) {
            log << strformat("Map %ld to %ld (parent %ld nchild %d)", i, map[i], parent[i], children[i].nChild)
                << std::endl;
        }

    /* Set parents, children, root */
    this->root = map[root];
    int64_t node;
    for (node = 0; node < maxnodes; node++) {
        int64_t njnode = map[node];
        if (njnode >= 0) {
            child[njnode].nChild = children[node].nChild;
            for (int64_t i = 0; i < children[node].nChild; i++) {
                assert(children[node].child[i] >= 0 && children[node].child[i] < maxnodes);
                child[njnode].child[i] = map[children[node].child[i]];
            }
            if (parent[node] >= 0) {
                this->parent[njnode] = map[parent[node]];
            }
        }
    }

    /* Make sure that parents/child relationships match */
    for (int64_t i = 0; i < this->maxnode; i++) {
        Children &c = child[i];
        for (int64_t j = 0; j < c.nChild; j++) {
            assert(c.child[j] >= 0 && c.child[j] < this->maxnode && this->parent[c.child[j]] == i);
        }
    }
    assert(this->parent[this->root] < 0);

    /* Compute profiles as balanced -- the NNI stage will recompute these
       profiles anyway
    */
    std::vector<ibool> traversal(this->maxnode, false);
    node = root;
    while ((node = traversePostorder(node, traversal, nullptr, root)) >= 0) {//not parallelize
        if (node >= (int64_t) nSeqs && node != root) {
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
    constexpr auto FP_FORMAT = sizeof(numeric_t) == sizeof(float) ? "%.5f" : "%.9f";

    if (nSeqs == 1 && unique.alnNext[unique.uniqueFirst[0]] >= 0) {
        /* Special case -- otherwise we end up with double parens */
        int64_t first = unique.uniqueFirst[0];
        assert(first >= 0 && first < (int64_t) unique.alnToUniq.size());
        out << "(";
        quotes(out, names[first], options.bQuote) << ":0.0";
        int64_t iName = unique.alnNext[first];
        while (iName >= 0) {
            assert(iName < (int64_t) unique.alnToUniq.size());
            out << ",";
            quotes(out, names[iName], options.bQuote) << ":0.0";
            iName = unique.alnNext[iName];
        }
        out << ");" << std::endl;
        return;
    }

    int64_t stackSize = 1;
    std::vector<std::pair<int64_t, int64_t>> stack(maxnodes);
    stack[0].first = root; //node
    stack[0].second = 0; //end

    while (stackSize > 0) {
        std::pair<int64_t, int64_t> &last = stack[stackSize - 1];
        stackSize--;
        /* Save last, as we are about to overwrite it */
        int64_t node = last.first;
        int64_t end = last.second;

        if (node < (int64_t) nSeqs) {
            if (child[parent[node]].child[0] != node) {
                out << ",";
            }
            int64_t first = unique.uniqueFirst[node];
            assert(first >= 0 && first < (int64_t) unique.alnToUniq.size());
            /* Print the name, or the subtree of duplicate names */
            if (unique.alnNext[first] == -1) {
                quotes(out, names[first], options.bQuote);
            } else {
                out << "(";
                quotes(out, names[first], options.bQuote) << ":0.0";
                int64_t iName = unique.alnNext[first];
                while (iName >= 0) {
                    assert(iName < (int64_t) unique.alnToUniq.size());
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
            } else if (bShowSupport) {
                out << strformat(")%.3f:", support[node]) << strformat(FP_FORMAT, branchlength[node]);
            } else {
                out << "):" << strformat(FP_FORMAT, branchlength[node]);
            }
        } else {
            if (node != root && child[parent[node]].child[0] != node) {
                out << ",";
            }
            out << "(";
            stackSize++;
            stack[stackSize - 1].first = node;
            stack[stackSize - 1].second = 1;
            Children &c = child[node];
            /* put children on in reverse order because we use the last one first */
            for (int i = c.nChild - 1; i >= 0; i--) {
                stackSize++;
                stack[stackSize - 1].first = c.child[i];
                stack[stackSize - 1].second = 0;
            }
        }
    }
    out << ";" << std::endl;
}

AbsNeighbourJoining(void)::fastNJ() {
    assert(nSeqs >= 1);
    if (nSeqs < 3) {
        root = maxnode++;
        child[root].nChild = (int) nSeqs;
        for (int64_t iNode = 0; iNode < (int64_t) nSeqs; iNode++) {
            parent[iNode] = root;
            child[root].child[iNode] = iNode;
        }
        if (nSeqs == 1) {
            branchlength[0] = 0;
        } else {
            assert (nSeqs == 2);
            Besthit hit;
            seqDist(profiles[0].codes, profiles[1].codes, hit);
            branchlength[0] = hit.dist / 2.0;
            branchlength[1] = hit.dist / 2.0;
        }
        return;
    }

    /* else 3 or more sequences */

    /* The visible set stores the best hit of each node (unless using top hits, in which case
       it is handled by the top hits routines) */
    std::vector<Besthit> visible;    /* Not used if doing top hits */
    std::vector<Besthit> besthitNew; /* All hits of new node -- not used if doing top-hits */

    /* The top-hits lists, with the key parameter m = length of each top-hit list */
    std::unique_ptr<TopHits> tophits;
    int64_t m = 0;            /* maximum length of a top-hits list */
    if (options.tophitsMult > 0) {
        m = (int64_t) (0.5 + options.tophitsMult * sqrt(nSeqs));
        if (m < 4 || 2 * m >= (int64_t) nSeqs) {
            m = 0;
            if (options.verbose > 1) {
                log << "Too few leaves, turning off top-hits" << std::endl;
            }
        } else {
            if (options.verbose > 2) {
                log << strformat("Top-hit-list size = %ld of %ld", m, nSeqs) << std::endl;
            }
        }
    }
    assert(!(options.slow && m > 0));

    /* Initialize top-hits or visible set */
    if (m > 0) {
        tophits = make_unique<TopHits>(options, maxnodes, m);
        setAllLeafTopHits(*tophits);
        resetTopVisible((int64_t) nSeqs, *tophits);
    } else if (!options.slow) {
        visible.resize(maxnodes);
        besthitNew.resize(maxnodes);
        for (int64_t iNode = 0; iNode < (int64_t) nSeqs; iNode++) {
            setBestHit(iNode, /*nActive*/nSeqs, visible[iNode], /*OUT IGNORED*/nullptr);
        }
    }

    /* Iterate over joins */
    int64_t nActiveOutProfileReset = nSeqs;
    for (int64_t nActive = nSeqs; nActive > 3; nActive--) {
        int64_t nJoinsDone = nSeqs - nActive;
        if (nJoinsDone > 0 && (nJoinsDone % 100) == 0) {
            progressReport.print("Joined %6ld of %6ld", nJoinsDone, (int64_t) (nSeqs - 3));
        }

        Besthit join;        /* the join to do */
        if (options.slow) {
            exhaustiveNJSearch(nActive,/*OUT*/join);
        } else if (m > 0) {
            topHitNJSearch(nActive, /*IN/OUT*/*tophits, /*OUT*/join);
        } else {
            fastNJSearch(nActive, /*IN/OUT*/visible, /*OUT*/join);
        }

        if (options.verbose > 2) {
            double penalty = options.constraintWeight * (double) joinConstraintPenalty(join.i, join.j);
            if (penalty > 0.001) {
                log << strformat("Constraint violation during neighbor-joining %ld %ld into %ld penalty %.3f",
                                 join.i, join.j, maxnode, penalty) << std::endl;

                for (int64_t iC = 0; iC < (int64_t) constraintSeqs.size(); iC++) {
                    int64_t local = joinConstraintPenaltyPiece(join.i, join.j, iC);
                    if (local > 0)
                        log << strformat("Constraint %ld piece %ld %ld/%ld %ld/%ld %ld/%ld", iC, local,
                                         profiles[join.i].nOn[iC],
                                         profiles[join.i].nOff[iC],
                                         profiles[join.j].nOn[iC],
                                         profiles[join.j].nOff[iC],
                                         outprofile.nOn[iC] - profiles[join.i].nOn[iC] - profiles[join.j].nOn[iC],
                                         outprofile.nOff[iC] - profiles[join.i].nOff[iC] -
                                         profiles[join.j].nOff[iC]) << std::endl;
                }
            }
        }

        /* because of the stale out-distance heuristic, make sure that these are up-to-date */
        setOutDistance(join.i, nActive);
        setOutDistance(join.j, nActive);
        /* Make sure weight is set and criterion is up to date */
        setDistCriterion(nActive, join);
        assert(nOutDistActive[join.i] == nActive);
        assert(nOutDistActive[join.j] == nActive);

        int64_t newnode = maxnode++;
        parent[join.i] = newnode;
        parent[join.j] = newnode;
        child[newnode].nChild = 2;
        child[newnode].child[0] = join.i < join.j ? join.i : join.j;
        child[newnode].child[1] = join.i > join.j ? join.i : join.j;

        double rawIJ = join.dist + diameter[join.i] + diameter[join.j];
        double distIJ = join.dist;

        double deltaDist = (outDistances[join.i] - outDistances[join.j]) / (double) (nActive - 2);
        branchlength[join.i] = (distIJ + deltaDist) / 2;
        branchlength[join.j] = (distIJ - deltaDist) / 2;

        double bionjWeight = 0.5;    /* IJ = bionjWeight*I + (1-bionjWeight)*J */
        double varIJ = rawIJ - varDiameter[join.i] - varDiameter[join.j];

        if (options.bionj && join.weight > 0.01 && varIJ > 0.001) {
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
            Besthit outI;
            Besthit outJ;
            profileDist(profiles[join.i], outprofile, outI);
            profileDist(profiles[join.j], outprofile, outJ);
            options.debug.outprofileOps += 2;

            double varIWeight = (nActive * outI.weight - selfweight[join.i] - join.weight);
            double varJWeight = (nActive * outJ.weight - selfweight[join.j] - join.weight);

            double varITop = outI.dist * outI.weight * nActive
                             - selfdist[join.i] * selfweight[join.i] - rawIJ * join.weight;
            double varJTop = outJ.dist * outJ.weight * nActive
                             - selfdist[join.j] * selfweight[join.j] - rawIJ * join.weight;

            double deltaProfileVarOut = (nActive - 2) * (varJTop / varJWeight - varITop / varIWeight);
            double deltaVarDiam = (nActive - 2) * (varDiameter[join.i] - varDiameter[join.j]);
            if (varJWeight > 0.01 && varIWeight > 0.01)
                bionjWeight = 0.5 + (deltaProfileVarOut + deltaVarDiam) / (2 * (nActive - 2) * varIJ);
            if (bionjWeight < 0) bionjWeight = 0;
            if (bionjWeight > 1) bionjWeight = 1;
            if (options.verbose > 2) {
                log << strformat("dVarO %f dVarDiam %f varIJ %f from dist %f weight %f (pos %ld) bionjWeight %f %f",
                                 deltaProfileVarOut, deltaVarDiam, varIJ, join.dist, join.weight, nPos, bionjWeight,
                                 1 - bionjWeight);
            }
            if (options.verbose > 3 && (newnode % 5) == 0) {
                /* Compare weight estimated from outprofiles from weight made by summing over other nodes */
                double deltaProfileVarTot = 0;
                for (int64_t iNode = 0; iNode < newnode; iNode++) {
                    if (parent[iNode] < 0) { /* excludes join.i, join.j */
                        Besthit di, dj;
                        profileDist(profiles[join.i], profiles[iNode], di);
                        profileDist(profiles[join.j], profiles[iNode], dj);
                        deltaProfileVarTot += dj.dist - di.dist;
                    }
                }
                double lambdaTot = 0.5 + (deltaProfileVarTot + deltaVarDiam) / (2 * (nActive - 2) * varIJ);
                if (lambdaTot < 0) {
                    lambdaTot = 0;
                }
                if (lambdaTot > 1) {
                    lambdaTot = 1;
                }
                if (fabs(bionjWeight - lambdaTot) > 0.01 || options.verbose > 4) {
                    fprintf(stderr, "deltaProfileVar actual %.6f estimated %.6f lambda actual %.3f estimated %.3f\n",
                            deltaProfileVarTot, deltaProfileVarOut, lambdaTot, bionjWeight);
                }
            }
        }
        if (options.verbose > 2) {
            log << strformat("Join\t%ld\t%ld\t%.6f\tlambda\t%.6f\tselfw\t%.3f\t%.3f\tnew\t%ld",
                             join.i < join.j ? join.i : join.j,
                             join.i < join.j ? join.j : join.i,
                             join.criterion, bionjWeight,
                             selfweight[join.i < join.j ? join.i : join.j],
                             selfweight[join.i < join.j ? join.j : join.i],
                             newnode) << std::endl;
        }

        diameter[newnode] = bionjWeight * (branchlength[join.i] + diameter[join.i])
                            + (1 - bionjWeight) * (branchlength[join.j] + diameter[join.j]);
        varDiameter[newnode] = bionjWeight * varDiameter[join.i]
                               + (1 - bionjWeight) * varDiameter[join.j]
                               + bionjWeight * (1 - bionjWeight) * varIJ;
        averageProfile(profiles[newnode], profiles[join.i], profiles[join.j],
                       options.bionj ? bionjWeight : /*noweight*/-1.0);

        /* Update out-distances and total diameters */
        int64_t changedActiveOutProfile = nActiveOutProfileReset - (nActive - 1);
        if (changedActiveOutProfile >= options.nResetOutProfile
            && changedActiveOutProfile >= options.fResetOutProfile * nActiveOutProfileReset) {
            /* Recompute the outprofile from scratch to avoid roundoff error */
            std::vector<Profile *> activeProfiles(nActive - 1);

            int64_t nSaved = 0;
            totdiam = 0;
            for (int64_t iNode = 0; iNode < maxnode; iNode++) {
                if (parent[iNode] < 0) {
                    assert(nSaved < nActive - 1);
                    activeProfiles[nSaved++] = &profiles[iNode];
                    totdiam += diameter[iNode];
                }
            }
            assert(nSaved == nActive - 1);
            if (options.verbose > 2) {
                log << strformat("Recomputing outprofile %ld %ld", nActiveOutProfileReset, nActive - 1) << std::endl;
            }
            outProfile(outprofile, activeProfiles, nSaved);
            nActiveOutProfileReset = nActive - 1;
        } else {
            updateOutProfile(/*OUT*/outprofile, profiles[join.i], profiles[join.j], profiles[newnode], nActive);
            totdiam += diameter[newnode] - diameter[join.i] - diameter[join.j];
        }

        /* Store self-dist for use in other computations */
        Besthit _selfdist;
        profileDist(profiles[newnode], profiles[newnode], _selfdist);
        selfdist[newnode] = _selfdist.dist;
        selfweight[newnode] = _selfdist.weight;

        /* Find the best hit of the joined node IJ */
        if (m > 0) {
            topHitJoin(newnode, nActive - 1, *tophits);
        } else {
            /* Not using top-hits, so we update all out-distances */
            for (int64_t iNode = 0; iNode < maxnode; iNode++) {
                if (parent[iNode] < 0) {
                    /* True nActive is now nActive-1 */
                    setOutDistance(iNode, nActive - 1);
                }
            }

            if (!visible.empty()) {
                setBestHit(newnode, nActive - 1, visible[newnode], /*OUT OPTIONAL*/besthitNew.data());
                if (options.verbose > 2) {
                    log << strformat("Visible %ld %ld %f %f",
                                     visible[newnode].i, visible[newnode].j,
                                     visible[newnode].dist, visible[newnode].criterion) << std::endl;
                }
                if (!besthitNew.empty()) {
                    /* Use distances to new node to update visible set entries that are non-optimal */
                    for (int64_t iNode = 0; iNode < maxnode; iNode++) {
                        if (parent[iNode] >= 0 || iNode == newnode) {
                            continue;
                        }
                        auto iOldVisible = visible[iNode].j;
                        assert(iOldVisible >= 0);
                        assert(visible[iNode].i == iNode);

                        /* Update the criterion; use nActive-1 because haven't decremented nActive yet */
                        if (parent[iOldVisible] < 0) {
                            setCriterion(nActive - 1, visible[iNode]);
                        }

                        if (parent[iOldVisible] >= 0
                            || besthitNew[iNode].criterion < visible[iNode].criterion) {
                            if (options.verbose > 3) {
                                log << strformat("Visible %ld reset from %ld to %ld (%f vs. %f)",
                                                 iNode, iOldVisible,
                                                 newnode, visible[iNode].criterion, besthitNew[iNode].criterion)
                                    << std::endl;
                            }
                            if (parent[iOldVisible] < 0) {
                                options.debug.nVisibleUpdate++;
                            }
                            visible[iNode].j = newnode;
                            visible[iNode].dist = besthitNew[iNode].dist;
                            visible[iNode].criterion = besthitNew[iNode].criterion;
                        }
                    } /* end loop over all nodes */
                } /* end if recording all hits of new node */
            } /* end if keeping a visible set */
        } /* end else (m==0) */
    } /* end loop over nActive */

    /* We no longer need the tophits, visible set, etc. */
    visible.clear();
    visible.reserve(0);
    besthitNew.clear();
    besthitNew.reserve(0);

    /* Add a root for the 3 remaining nodes */
    int64_t top[3];
    int64_t nTop = 0;
    for (int64_t iNode = 0; iNode < maxnode; iNode++) {
        if (parent[iNode] < 0) {
            assert(nTop <= 2);
            top[nTop++] = iNode;
        }
    }
    assert(nTop == 3);

    root = maxnode++;
    child[root].nChild = 3;
    for (nTop = 0; nTop < 3; nTop++) {
        parent[top[nTop]] = root;
        child[root].child[nTop] = top[nTop];
    }

    Besthit dist01, dist02, dist12;
    profileDist(profiles[top[0]], profiles[top[1]], dist01);
    profileDist(profiles[top[0]], profiles[top[2]], dist02);
    profileDist(profiles[top[1]], profiles[top[2]], dist12);

    double d01 = dist01.dist - diameter[top[0]] - diameter[top[1]];
    double d02 = dist02.dist - diameter[top[0]] - diameter[top[2]];
    double d12 = dist12.dist - diameter[top[1]] - diameter[top[2]];
    branchlength[top[0]] = (d01 + d02 - d12) / 2;
    branchlength[top[1]] = (d01 + d12 - d02) / 2;
    branchlength[top[2]] = (d02 + d12 - d01) / 2;

    /* Check how accurate the outprofile is */
    if (options.verbose > 2) {
        std::vector<Profile *> tmp = {&profiles[top[0]], &profiles[top[1]], &profiles[top[2]]};
        Profile out(nPos, constraintSeqs.size());
        /* Use swap to avoid a deep copy of Profile*/
        outProfile(out, tmp, 3);
        double freqerror = 0;
        double weighterror = 0;
        for (int64_t i = 0; i < nPos; i++) {
            weighterror += fabs(out.weights[i] - outprofile.weights[i]);
            for (int k = 0; k < options.nCodes; k++) {
                freqerror += fabs(out.vectors[nCodeSize * i + k] - outprofile.vectors[nCodeSize * i + k]);
            }
        }

        log << strformat("Roundoff error in outprofile@end: WeightError %f FreqError %f", weighterror, freqerror)
            << std::endl;
    }
}

AbsNeighbourJoining(int64_t)::traverseReliabilityNJ(int64_t node, const std::vector<int64_t> &col,
                                                    std::unique_ptr<Profile> upProfiles[],
                                                    std::vector<ibool> &traversal) {
    const bool parallel = omp_in_parallel();
    int64_t branchRoot = node;
    int64_t iNodesDone = 0;
    while ((node = traversePostorder(node,/*IN/OUT*/traversal, /*pUp*/nullptr, branchRoot)) >= 0) {//level-1
        if (node < (int64_t) nSeqs || node == root)
            continue; /* nothing to do for leaves or root */

        if (iNodesDone > 0 && (iNodesDone % 100) == 0 && !parallel)
            progressReport.print("Local bootstrap for %6ld of %6ld internal splits", iNodesDone, (int64_t) (nSeqs - 3));
        iNodesDone++;

        Profile *profiles4[4];
        int64_t nodeABCD[4];
        setupABCD(node, /*OUT*/profiles4, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, /*useML*/false);

        support[node] = splitSupport(*profiles4[0], *profiles4[1], *profiles4[2], *profiles4[3], col);
        /* no longer needed */
        upProfiles[nodeABCD[0]].reset();
        upProfiles[nodeABCD[1]].reset();
        upProfiles[nodeABCD[2]].reset();
    }
    return iNodesDone;
}

AbsNeighbourJoining(void)::reliabilityNJ() {
    /* For each non-root node N, with children A,B, parent P, sibling C, and grandparent G,
       we test the reliability of the split (A,B) versus rest by comparing the profiles
       of A, B, C, and the "up-profile" of P.

       Each node's upProfile is the average of its sibling's (down)-profile + its parent's up-profile
       (If node's parent is the root, then there are two siblings and we don't need an up-profile)

       To save memory, we do depth-first-search down from the root, and we only keep
       up-profiles for nodes in the active path.
    */
    if (nSeqs <= 3 || options.nBootstrap <= 0) {
        return;            /* nothing to do */
    }
    std::vector<int64_t> col;

    resampleColumns(col);

    std::vector<ibool> traversal(maxnodes, false);

    if (options.threads > 1 && options.threadsLevel > 0) {
        std::vector<std::vector<int64_t>> chunks;
        treePartition(chunks, traversal, 0);
        int64_t done = 0;

        #pragma omp parallel
        {
            std::vector<std::unique_ptr<Profile>> upProfiles2(maxnodes);
            int64_t done2;

            #pragma omp for schedule(static, 1)
            for (int i = 0; i < (int) chunks.size(); i++) {
                for (int64_t subtreeRoot: chunks[i]) {
                    done2 = traverseReliabilityNJ(subtreeRoot, col, upProfiles2.data(), traversal);
                    if (options.verbose > 0) {
                        #pragma omp critical
                        {
                            done += done2;
                            progressReport.print("Local bootstrap for %6ld of %6ld internal splits", done,
                                                 (int64_t) (nSeqs - 3));
                        }
                    }
                }
            }
        }
    }

    std::vector<std::unique_ptr<Profile>> upProfiles(maxnodes);
    traverseReliabilityNJ(root, col, upProfiles.data(), traversal);
}

AbsNeighbourJoining(void)::
readTreeAddChild(int64_t iparent, int64_t ichild, std::vector<int64_t> &parents, std::vector<Children> &children) {
    assert(iparent >= 0);
    assert(ichild >= 0);
    assert(parents[ichild] < 0);
    assert(children[iparent].nChild < 3);
    parents[ichild] = iparent;
    children[iparent].child[children[iparent].nChild++] = ichild;
}

AbsNeighbourJoining(void)::readTreeMaybeAddLeaf(int64_t iparent, std::string &name, HashTable &hashnames,
                                                Uniquify &unique,
                                                std::vector<int64_t> &parents, std::vector<Children> &children) {
    auto hi = hashnames.find(name);
    if (hi == nullptr) {
        readTreeError("not recognized as a sequence name", name);
    }

    int64_t iSeqNonunique = *hi;
    assert(iSeqNonunique >= 0 && iSeqNonunique < (int64_t) unique.alnToUniq.size());
    int64_t iSeqUnique = unique.alnToUniq[iSeqNonunique];
    assert(iSeqUnique >= 0 && iSeqUnique < (int64_t) unique.uniqueSeq.size());
    /* Either record this leaves' parent (if it is -1) or ignore this leaf (if already seen) */
    if (parents[iSeqUnique] < 0) {
        readTreeAddChild(iparent, iSeqUnique, /*IN/OUT*/parents, /*IN/OUT*/children);
        if (options.verbose > 5) {
            log << strformat("Found leaf uniq%ld name %s child of %ld", iSeqUnique, name.c_str(), iparent) << std::endl;
        }
    } else {
        if (options.verbose > 5) {
            log << strformat("Skipped redundant leaf uniq%ld name %s", iSeqUnique, name.c_str()) << std::endl;
        }
    }
}

AbsNeighbourJoining(void)::readTreeRemove(std::vector<int64_t> &parents, std::vector<Children> &children,
                                          int64_t node) {
    if (options.verbose > 5) {
        log << strformat("Removing node %ld parent %ld", node, parents[node]) << std::endl;
    }
    assert(parents[node] >= 0);
    int64_t iparent = parents[node];
    parents[node] = -1;
    Children &pc = children[iparent];
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
                log << strformat("Repointing parent %ld to child %ld", iparent, nc.child[j]) << std::endl;
            }
            pc.child[pc.nChild++] = nc.child[j];
            parents[nc.child[j]] = iparent;
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
            if (buf.empty()) {
                buf += (char) c;
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
            buf += (char) c;
        }
    }
    return !buf.empty();
}

AbsNeighbourJoining(int64_t)::
traversePostorder(int64_t node, std::vector<ibool> &traversal, bool *pUp, int64_t branchRoot) {
    if (pUp != nullptr) {
        *pUp = false;
    }
    while (true) {
        assert(node >= 0);

        /* move to a child if possible */
        bool found = false;
        int iChild;
        for (iChild = 0; iChild < child[node].nChild; iChild++) {
            int64_t childnode = child[node].child[iChild];
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
        if (node == branchRoot) {
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

AbsNeighbourJoining(typename veryfasttree::NeighbourJoining<Precision, Operations>::Profile *)::getUpProfile(
        std::unique_ptr<Profile> upProfiles[], int64_t outnode, bool useML) {
    assert(outnode != root && outnode >= (int64_t) nSeqs); /* not for root or leaves */
    if (upProfiles[outnode]) {
        return upProfiles[outnode].get();
    }

    std::vector<int64_t> path;
    pathToRoot(outnode, path);

    /* depth-1 is root */
    for (int64_t i = path.size() - 2; i >= 0; i--) {
        int64_t node = path[i];

        if (!upProfiles[node]) {
            /* Note -- SetupABCD may call GetUpProfile, but it should do it farther
               up in the path to the root
            */
            Profile *profiles4[4];
            int64_t nodeABCD[4];
            setupABCD(node, profiles4, upProfiles, nodeABCD, useML);
            upProfiles[node] = make_unique<Profile>(nPos, constraintSeqs.size());
            if (useML) {
                /* If node is a child of root, then the 4th profile is of the 2nd root-sibling of node
                   Otherwise, the 4th profile is the up-profile of the parent of node, and that
                   is the branch-length we need
                 */
                double lenC = branchlength[nodeABCD[2]];
                double lenD = branchlength[nodeABCD[3]];
                if (options.verbose > 3) {
                    log << strformat("Computing UpProfile for node %ld with lenC %.4f lenD %.4f pair-loglk %.3f",
                                     node, lenC, lenD,
                                     pairLogLk(*profiles4[2], *profiles4[3], lenC + lenD, /*site_lk*/ nullptr))
                        << std::endl;
                    printNJInternal(log, /*useLen*/true);
                }
                posteriorProfile(*upProfiles[node], /*C*/*profiles4[2], /*D*/*profiles4[3], lenC, lenD);
            } else {
                Profile *profiles4CDAB[4] = {profiles4[2], profiles4[3], profiles4[0], profiles4[1]};
                double weight = quartetWeight(profiles4CDAB);
                if (options.verbose > 3) {
                    log << strformat(
                            "Compute upprofile of %ld from %ld and parents (vs. children %ld %ld) with weight %.3f",
                            node, nodeABCD[2], nodeABCD[0], nodeABCD[1], weight) << std::endl;
                }
                averageProfile(*upProfiles[node], *profiles4[2], *profiles4[3], weight);
            }
        }
    }
    assert(upProfiles[outnode]);
    return upProfiles[outnode].get();
}

AbsNeighbourJoining(void)::recomputeProfile(std::unique_ptr<Profile> upProfiles[], int64_t node, int64_t useML) {
    if (node < (int64_t) nSeqs || node == root) {
        return;            /* no profile to compute */
    }
    assert(child[node].nChild == 2);

    Profile *profiles4[4] = {nullptr, nullptr, nullptr, nullptr};
    double weight = 0.5;
    if (useML || !options.bionj) {
        profiles4[0] = &profiles[child[node].child[0]];
        profiles4[1] = &profiles[child[node].child[1]];
    } else {
        int64_t nodeABCD[4];
        setupABCD(node, profiles4, upProfiles, nodeABCD, useML);
        weight = quartetWeight(profiles4);
    }
    if (options.verbose > 3) {
        if (useML) {
            log << strformat("Recompute %ld from %ld %ld lengths %.4f %.4f",
                             node,
                             child[node].child[0],
                             child[node].child[1],
                             branchlength[child[node].child[0]],
                             branchlength[child[node].child[1]]) << std::endl;
        } else {
            log << strformat("Recompute %ld from %ld %ld weight %.3f",
                             node, child[node].child[0], child[node].child[1], weight) << std::endl;
        }
    }
    if (useML) {
        posteriorProfile(profiles[node], *profiles4[0], *profiles4[1],
                         branchlength[child[node].child[0]],
                         branchlength[child[node].child[1]]);
    } else {
        averageProfile(profiles[node], *profiles4[0], *profiles4[1], weight);
    }
}

AbsNeighbourJoining(inline void)::traverseRecomputeProfiles(int64_t node, std::vector<ibool> &traversal,
                                                            DistanceMatrix <Precision, op_t::ALIGNMENT> &dmat) {
    int64_t branchRoot = node;
    while ((node = traversePostorder(node, traversal, nullptr, branchRoot)) >= 0) {//level-1
        if (child[node].nChild == 2) {
            int64_t *children = child[node].child;
            averageProfile(profiles[node], profiles[children[0]], profiles[children[1]],  /*unweighted*/-1.0, dmat);
        }
    }
}

AbsNeighbourJoining(void)::recomputeProfiles(DistanceMatrix <Precision, op_t::ALIGNMENT> &dmat) {
    std::vector<ibool> traversal(maxnodes, false);

    if (options.threads > 1 && options.threadsLevel > 0) {
        std::vector<std::vector<int64_t>> chunks;
        treePartition(chunks, traversal, 0);

        #pragma omp parallel
        {
            #pragma omp for schedule(static, 1)
            for (int i = 0; i < (int) chunks.size(); i++) {
                for (int64_t subtreeRoot: chunks[i]) {
                    traverseRecomputeProfiles(subtreeRoot, traversal, dmat);
                }
            }
        }
    }

    traverseRecomputeProfiles(root, traversal, dmat);
}


AbsNeighbourJoining(inline void)::traverseRecomputeMLProfiles(int64_t node, std::vector<ibool> &traversal) {
    int64_t branchRoot = node;
    while ((node = traversePostorder(node, traversal, nullptr, branchRoot)) >= 0) {//level-1
        if (child[node].nChild == 2) {
            int64_t *children = child[node].child;
            posteriorProfile(profiles[node], profiles[children[0]], profiles[children[1]],
                             branchlength[children[0]], branchlength[children[1]]);
        }
    }
}

AbsNeighbourJoining(void)::recomputeMLProfiles() {
    std::vector<ibool> traversal(maxnodes, false);

    if (options.threads > 1 && options.threadsLevel > 0) {
        std::vector<std::vector<int64_t>> chunks;
        treePartition(chunks, traversal, 0);

        #pragma omp parallel
        {
            #pragma omp for schedule(static, 1)
            for (int i = 0; i < (int) chunks.size(); i++) {
                for (int64_t subtreeRoot: chunks[i]) {
                    traverseRecomputeMLProfiles(subtreeRoot, traversal);
                }
            }
        }
    }

    traverseRecomputeMLProfiles(root, traversal);

}

/* The BIONJ-like formula for the weight of A when building a profile for AB is
     1/2 + (avgD(B,CD) - avgD(A,CD))/(2*d(A,B))
*/
AbsNeighbourJoining(double)::quartetWeight(Profile *profiles4[4]) {
    if (!options.bionj) {
        return (-1.0); /* even weighting */
    }
    double d[6];
    correctedPairDistances(profiles4, 4, d);
    if (d[qAB] < 0.01) {
        return -1.0;
    }
    double weight = 0.5 + ((d[qBC] + d[qBD]) - (d[qAC] + d[qAD])) / (4 * d[qAB]);
    if (weight < 0) {
        weight = 0;
    }
    if (weight > 1) {
        weight = 1;
    }
    return weight;
}

AbsNeighbourJoining(void)::pathToRoot(int64_t node, std::vector<int64_t> &path) {
    int64_t ancestor = node;
    while (ancestor >= 0) {
        path.push_back(ancestor);
        ancestor = parent[ancestor];
    }
}

AbsNeighbourJoining(void)::setBestHit(int64_t node, int64_t nActive, Besthit &bestjoin, Besthit allhits[],
                                      bool shared) {
    assert(parent[node] < 0);

    if (shared) {
        bestjoin.i = node;
        bestjoin.j = -1;
        bestjoin.dist = (numeric_t) 1e20;
        bestjoin.criterion = (numeric_t) 1e20;
        #pragma omp barrier
        #pragma omp for schedule(dynamic)
        for (int64_t j = 0; j < maxnode; j++) {
            Besthit &sv = allhits[j];
            sv.i = node;
            sv.j = j;
            if (parent[j] >= 0) {
                sv.i = -1;        /* illegal/empty join */
                sv.weight = 0.0;
                sv.criterion = sv.dist = (numeric_t) 1e20;
                continue;
            }
            /* Note that we compute self-distances (allow j==node) because the top-hit heuristic
               expects self to be within its top hits, but we exclude those from the bestjoin
               that we return...
            */
            setDistCriterion(nActive, sv);
        }
    } else {
        bestjoin.i = node;
        bestjoin.j = -1;
        bestjoin.dist = (numeric_t) 1e20;
        bestjoin.criterion = (numeric_t) 1e20;

        Besthit tmp;

        #pragma omp parallel
        {
            Besthit bestjoin2 = bestjoin;
            #pragma omp for schedule(dynamic)
            for (int64_t j = 0; j < maxnode; j++) {
                Besthit &sv = allhits != nullptr ? allhits[j] : tmp;
                sv.i = node;
                sv.j = j;
                if (parent[j] >= 0) {
                    sv.i = -1;        /* illegal/empty join */
                    sv.weight = 0.0;
                    sv.criterion = sv.dist = (numeric_t) 1e20;
                    continue;
                }
                /* Note that we compute self-distances (allow j==node) because the top-hit heuristic
                   expects self to be within its top hits, but we exclude those from the bestjoin
                   that we return...
                */
                setDistCriterion(nActive, sv);
                /*
                 * Local comparisons in each thread to remove non deterministic results
                 */
                if (sv.criterion < bestjoin2.criterion && node != j) {
                    bestjoin2 = sv;
                }
            }

            #pragma omp critical
            {
                if (bestjoin2.criterion < bestjoin.criterion) {
                    bestjoin = bestjoin2;
                }
            }
        }
    }

    if (options.verbose > 5) {
        log << strformat("SetBestHit %ld %ld %f %f", bestjoin.i, bestjoin.j, bestjoin.dist, bestjoin.criterion)
            << std::endl;
    }
}

AbsNeighbourJoining(void)::exhaustiveNJSearch(int64_t nActive, Besthit &join) {
    join.i = -1;
    join.j = -1;
    join.weight = 0;
    join.dist = (numeric_t) 1e20;
    join.criterion = (numeric_t) 1e20;

    #pragma omp parallel
    {
        Besthit join2 = join;
        #pragma omp  for schedule(dynamic)
        for (int64_t i = 0; i < maxnode - 1; i++) {
            if (parent[i] < 0) {
                for (int64_t j = i + 1; j < maxnode; j++) {
                    if (parent[j] < 0) {
                        Besthit hit;
                        hit.i = i;
                        hit.j = j;
                        setDistCriterion(nActive, hit);
                        if (hit.criterion < join2.criterion) {
                            join2 = hit;
                        }
                    }
                }
            }
        }
        #pragma omp critical
        {
            if (join2.criterion < join.criterion) {
                join = join2;
            }
        };
    };

    assert (join.i >= 0 && join.j >= 0);
    std::cerr << join.i << " " << join.j << std::endl;
}

AbsNeighbourJoining(void)::fastNJSearch(int64_t nActive, std::vector<Besthit> &besthits, Besthit &join) {
    join.i = -1;
    join.j = -1;
    join.dist = (numeric_t) 1e20;
    join.weight = 0;
    join.criterion = (numeric_t) 1e20;

    for (int64_t iNode = 0; iNode < maxnode; iNode++) {
        int64_t jNode = besthits[iNode].j;
        if (parent[iNode] < 0 && parent[jNode] < 0) { /* both i and j still active */
            /* recompute criterion to reflect the current out-distances */
            setCriterion(nActive, besthits[iNode]);
            if (besthits[iNode].criterion < join.criterion) {
                join = besthits[iNode];
            }
        }
    }

    if (!options.fastest) {
        int changed;
        do {
            changed = 0;
            assert(join.i >= 0 && join.j >= 0);
            setBestHit(join.i, nActive, besthits[join.i], nullptr);
            if (besthits[join.i].j != join.j) {
                changed = 1;
                if (options.verbose > 2) {
                    log << strformat("BetterI\t%ld\t%ld\t%ld\t%ld\t%f\t%f",
                                     join.i, join.j, besthits[join.i].i, besthits[join.i].j,
                                     join.criterion, besthits[join.i].criterion) << std::endl;
                }
            }

            /* Save the best hit either way, because the out-distance has probably changed
               since we started the computation. */
            join.j = besthits[join.i].j;
            join.weight = besthits[join.i].weight;
            join.dist = besthits[join.i].dist;
            join.criterion = besthits[join.i].criterion;

            setBestHit(join.j, nActive, besthits[join.j], nullptr);
            if (besthits[join.j].j != join.i) {
                changed = 1;
                if (options.verbose > 2) {
                    log << strformat("BetterJ\t%ld\t%ld\t%ld\t%ld\t%f\t%f",
                                     join.i, join.j, besthits[join.j].i, besthits[join.j].j,
                                     join.criterion, besthits[join.j].criterion) << std::endl;
                }
                join.i = besthits[join.j].j;
                join.weight = besthits[join.j].weight;
                join.dist = besthits[join.j].dist;
                join.criterion = besthits[join.j].criterion;
            }
            if (changed) {
                options.debug.nHillBetter++;
            }
        } while (changed);
    }
}

AbsNeighbourJoining(void)::setAllLeafTopHits(TopHits &tophits_g) {
    double close = options.tophitsClose;
    if (close < 0) {
        if (options.fastest && nSeqs >= 50000) {
            close = 0.99;
        } else {
            double logN = std::log((double) nSeqs) / std::log(2.0);
            close = logN / (logN + 2.0);
        }
    }
    /* Sort the potential seeds, by a combination of nGaps and NJ->outDistances
       We don't store nGaps so we need to compute that
    */
    std::vector<int64_t> nGaps(nSeqs);

    for (int64_t iNode = 0; iNode < (int64_t) nSeqs; iNode++) {
        nGaps[iNode] = (int64_t) (0.5 + nPos - selfweight[iNode]);
    }

    std::vector<int64_t> seeds(nSeqs);
    for (int64_t iNode = 0; iNode < (int64_t) nSeqs; iNode++) {
        seeds[iNode] = iNode;
    }

    psort(seeds.begin(), seeds.end(), CompareSeeds(outDistances, nGaps));

    /* For each seed, save its top 2*m hits and then look for close neighbors */
    assert(2 * tophits_g.m <= (int64_t) nSeqs);

    int64_t nHasTopHits = 0;
    int64_t count = 0;
    std::vector<ibool> visited(nSeqs, false);

    if (options.deterministic && options.threads > 1) {
        std::vector<Besthit> besthitsSeed(nSeqs);
        TopHits &tophits = tophits_g;
        Besthit bestjoin;

        #pragma omp parallel firstprivate(nHasTopHits)
        {
            bool showlog = false;
            for (int64_t iSeed = 0; iSeed < (int64_t) nSeqs; iSeed++) {
                int64_t seed = seeds[iSeed];
                if (options.verbose > 0) {
                    if (iSeed % options.threads == omp_get_thread_num()) {
                        count += nHasTopHits;
                        nHasTopHits = 0;
                    }
                    if (iSeed % 100 == 0) {
                        showlog = true;
                    }
                }
                if (visited[seed]) {
                    continue;
                }
                setBestHit(seed, nSeqs, bestjoin, besthitsSeed.data(), true);
                #pragma omp master
                {
                    if (showlog) {
                        showlog = false;
                        progressReport.print("Top hits for %6ld of %6ld seqs (at seed %6ld)",
                                             count + 1, (int64_t) nSeqs, iSeed);
                    }
                    nHasTopHits++;
                    psort(besthitsSeed.begin(), besthitsSeed.end(), CompareHitsByCriterion(), true);
                    sortSaveBestHits(seed, besthitsSeed, nSeqs, tophits_g.m, tophits_g, false);
                };
                #pragma omp barrier
                visited[seed] = true;


                /* find "close" neighbors and compute their top hits */
                double neardist = besthitsSeed[2 * tophits.m - 1].dist * close;
                /* must have at least average weight, rem higher is better
                   and allow a bit more than average, e.g. if we are looking for within 30% away,
                   20% more gaps than usual seems OK
                   Alternatively, have a coverage requirement in case neighbor is short
                   If fastest, consider the top q/2 hits to be close neighbors, regardless
                */
                double nearweight = 0;
                for (int64_t iClose = 0; iClose < 2 * tophits.m; iClose++) {
                    nearweight += besthitsSeed[iClose].weight;
                }
                nearweight = nearweight / (2.0 * tophits.m); /* average */
                nearweight *= (1.0 - 2.0 * neardist / 3.0);
                double nearcover = 1.0 - neardist / 2.0;

                if (options.verbose > 2) {
                    log << strformat("Distance limit for close neighbors %f weight %f ungapped %ld",
                                     neardist, nearweight, nPos - nGaps[seed]) << std::endl;
                }

                #pragma omp  for schedule(dynamic)
                for (int64_t iClose = 0; iClose < tophits.m; iClose++) {
                    Besthit &closehit = besthitsSeed[iClose];
                    auto closeNode = closehit.j;
                    if (visited[closeNode]) {
                        continue;
                    }

                    /* If within close-distance, or identical, use as close neighbor */
                    bool isClose = closehit.dist <= neardist && (closehit.weight >= nearweight ||
                                                                 closehit.weight >=
                                                                 (nPos - nGaps[closeNode]) * nearcover);
                    bool identical = closehit.dist < 1e-6
                                     && fabs(closehit.weight - (nPos - nGaps[seed])) < 1e-5
                                     && fabs(closehit.weight - (nPos - nGaps[closeNode])) < 1e-5;
                    if (options.useTopHits2nd && iClose < tophits.q && (isClose || identical)) {
                        nHasTopHits++;
                        options.debug.nClose2Used++;
                        auto nUse = std::min(tophits.q * options.tophits2Safety, 2 * tophits.m);
                        std::vector<Besthit> besthitsClose(nUse);
                        transferBestHits(nSeqs, closeNode, besthitsSeed, nUse, besthitsClose.data(), true);
                        visited[closeNode] = true;
                        sortSaveBestHits(closeNode, besthitsClose, nUse, tophits.q, tophits);
                        tophits.topHitsLists[closeNode].hitSource = seed;
                    } else if (isClose || identical || (options.fastest && iClose < (tophits.q + 1) / 2)) {
                        nHasTopHits++;
                        options.debug.nCloseUsed++;
                        if (options.verbose > 2) {
                            log << strformat("Near neighbor %ld (rank %ld weight %f ungapped %ld %ld)",
                                             closeNode, iClose, besthitsSeed[iClose].weight,
                                             nPos - nGaps[seed],
                                             nPos - nGaps[closeNode]) << std::endl;
                        }

                        /* compute top 2*m hits */
                        std::vector<Besthit> besthitsNeighbor(2 * tophits.m);
                        transferBestHits(nSeqs, closeNode, besthitsSeed, 2 * tophits.m, besthitsNeighbor.data(),
                                         true);
                        visited[closeNode] = true;
                        sortSaveBestHits(closeNode, besthitsNeighbor, 2 * tophits.m, tophits.m, tophits);
                    }
                } /* end loop over close candidates */
            } /* end loop over seeds */
        }
    } else {
        std::vector<TopHits> threadTophits(0);
        {
            TopHits tophits_ref = tophits_g;
            tophits_ref.visible.resize(0);
            tophits_ref.topvisible.resize(0);
            threadTophits.resize(options.threads - 1, tophits_ref);
        }

        #pragma omp parallel firstprivate(nHasTopHits) if(!options.deterministic)
        {
            #pragma omp for schedule(static, (nSeqs / options.threads) + 1)
            for (int64_t iSeed = 0; iSeed < (int64_t) nSeqs; iSeed++) {
                int64_t seed = seeds[iSeed];
                if (options.threads == 1) {
                    if (iSeed > 0 && (iSeed % 100) == 0) {
                        progressReport.print("Top hits for %6ld of %6ld seqs (at seed %6ld)",
                                             nHasTopHits, (int64_t) nSeqs, iSeed);
                    }
                } else if (options.verbose > 0 && nHasTopHits > 100) {
                    #pragma omp critical
                    {
                        count += nHasTopHits;
                        progressReport.print("Top hits for %6ld of %6ld seqs (at seed %6ld)",
                                             count + 1, (int64_t) nSeqs, iSeed);
                        nHasTopHits = 0;
                    }
                }
                if (visited[seed]) {
                    if (options.verbose > 2) {
                        log << strformat("Skipping seed %ld", seed) << std::endl;
                    }
                    continue;
                }
                visited[seed] = true;
                TopHits &tophits = omp_get_thread_num() > 0 ? threadTophits[omp_get_thread_num() - 1] : tophits_g;

                std::vector<Besthit> besthitsSeed(nSeqs);
                std::vector<Besthit> besthitsNeighbor(2 * tophits.m);
                Besthit bestjoin;

                if (options.verbose > 2) {
                    log << strformat("Trying seed %ld", seed) << std::endl;
                }
                setBestHit(seed, nSeqs, bestjoin, besthitsSeed.data());

                /* sort & save top hits of self. besthitsSeed is now sorted. */
                sortSaveBestHits(seed, besthitsSeed, nSeqs, tophits.m, tophits);
                nHasTopHits++;

                /* find "close" neighbors and compute their top hits */
                double neardist = besthitsSeed[2 * tophits.m - 1].dist * close;
                /* must have at least average weight, rem higher is better
                   and allow a bit more than average, e.g. if we are looking for within 30% away,
                   20% more gaps than usual seems OK
                   Alternatively, have a coverage requirement in case neighbor is short
                   If fastest, consider the top q/2 hits to be close neighbors, regardless
                */
                double nearweight = 0;
                for (int64_t iClose = 0; iClose < 2 * tophits.m; iClose++) {
                    nearweight += besthitsSeed[iClose].weight;
                }
                nearweight = nearweight / (2.0 * tophits.m); /* average */
                nearweight *= (1.0 - 2.0 * neardist / 3.0);
                double nearcover = 1.0 - neardist / 2.0;

                if (options.verbose > 2) {
                    log << strformat("Distance limit for close neighbors %f weight %f ungapped %ld",
                                     neardist, nearweight, nPos - nGaps[seed]) << std::endl;
                }
                for (int64_t iClose = 0; iClose < tophits.m; iClose++) {
                    Besthit &closehit = besthitsSeed[iClose];
                    auto closeNode = closehit.j;
                    if (visited[closeNode]) {
                        continue;
                    }

                    /* If within close-distance, or identical, use as close neighbor */
                    bool isClose = closehit.dist <= neardist && (closehit.weight >= nearweight ||
                                                                 closehit.weight >=
                                                                 (nPos - nGaps[closeNode]) * nearcover);
                    bool identical = closehit.dist < 1e-6
                                     && fabs(closehit.weight - (nPos - nGaps[seed])) < 1e-5
                                     && fabs(closehit.weight - (nPos - nGaps[closeNode])) < 1e-5;
                    if (options.useTopHits2nd && iClose < tophits.q && (isClose || identical)) {
                        nHasTopHits++;
                        options.debug.nClose2Used++;
                        auto nUse = std::min(tophits.q * options.tophits2Safety, 2 * tophits.m);
                        std::vector<Besthit> besthitsClose(nUse);
                        transferBestHits(nSeqs, closeNode, besthitsSeed, nUse, besthitsClose.data(), true);
                        visited[closeNode] = true;
                        sortSaveBestHits(closeNode, besthitsClose, nUse, tophits.q, tophits);
                        tophits.topHitsLists[closeNode].hitSource = seed;
                    } else if (isClose || identical || (options.fastest && iClose < (tophits.q + 1) / 2)) {
                        nHasTopHits++;
                        options.debug.nCloseUsed++;
                        if (options.verbose > 2) {
                            log << strformat("Near neighbor %ld (rank %ld weight %f ungapped %ld %ld)",
                                             closeNode, iClose, besthitsSeed[iClose].weight,
                                             nPos - nGaps[seed],
                                             nPos - nGaps[closeNode]) << std::endl;
                        }

                        /* compute top 2*m hits */
                        transferBestHits(nSeqs, closeNode, besthitsSeed, 2 * tophits.m, besthitsNeighbor.data(),
                                         true);
                        visited[closeNode] = true;
                        sortSaveBestHits(closeNode, besthitsNeighbor, 2 * tophits.m, tophits.m, tophits);

                        /* And then try for a second level of transfer. We assume we
                           are in a good area, because of the 1st
                           level of transfer, and in a small neighborhood, because q is
                           small (32 for 1 million sequences), so we do not make any close checks.
                         */
                        for (int64_t iClose2 = 0; iClose2 < tophits.q && iClose2 < 2 * tophits.m; iClose2++) {
                            int64_t closeNode2 = besthitsNeighbor[iClose2].j;
                            assert(closeNode2 >= 0);
                            if (!visited[closeNode2]) {
                                options.debug.nClose2Used++;
                                nHasTopHits++;
                                auto nUse = std::min(tophits.q * options.tophits2Safety, 2 * tophits.m);
                                std::vector<Besthit> besthitsClose2(nUse);
                                transferBestHits(nSeqs, closeNode2, besthitsNeighbor, nUse, besthitsClose2.data(),
                                                 true);
                                visited[closeNode2] = true;
                                sortSaveBestHits(closeNode2, besthitsClose2, nUse, tophits.q, tophits);
                                tophits.topHitsLists[closeNode2].hitSource = closeNode;
                            } /* end if should do 2nd-level transfer */
                        }
                    }
                } /* end loop over close candidates */
            } /* end loop over seeds */

            // merge parallel tophits
            if (options.threads > 1) {
                #pragma omp  for schedule(dynamic)
                for (int64_t iSeed = 0; iSeed < (int64_t) nSeqs; iSeed++) {
                    for (int i = 0; i < options.threads - 1; i++) {
                        if (!tophits_g.topHitsLists[iSeed].hits.empty()) {
                            break;
                        }
                        if (!threadTophits[i].topHitsLists[iSeed].hits.empty()) {
                            std::swap(tophits_g.topHitsLists[iSeed].hits, threadTophits[i].topHitsLists[iSeed].hits);
                            break;
                        }
                    }
                }
            }
        }
    }

    TopHits &tophits = tophits_g;

    for (int64_t iNode = 0; iNode < (int64_t) nSeqs; iNode++) {
        TopHitsList &l = tophits.topHitsLists[iNode];
        assert(!l.hits.empty());
        assert(l.hits[0].j >= 0);
        assert(l.hits[0].j < (int64_t) nSeqs);
        assert(l.hits[0].j != iNode);
        tophits.visible[iNode] = l.hits[0];
    }

    if (options.verbose >= 2 && options.threads == 1) {
        log << strformat("#Close neighbors among leaves: 1st-level %ld 2nd-level %ld seeds %ld",
                         options.debug.nCloseUsed, options.debug.nClose2Used,
                         nSeqs - options.debug.nCloseUsed - options.debug.nClose2Used) << std::endl;
    }

    /* Now add a "checking phase" where we ensure that the q or 2*sqrt(m) hits
       of i are represented in j (if they should be)
     */
    int64_t lReplace = 0;
    int64_t nCheck = tophits.q > 0 ? tophits.q : (int64_t) (0.5 + 2.0 * sqrt(tophits.m));
    for (int64_t iNode = 0; iNode < (int64_t) nSeqs; iNode++) {
        if ((iNode % 100) == 0) {
            progressReport.print("Checking top hits for %6ld of %6ld seqs", iNode + 1, (int64_t) nSeqs);
        }
        TopHitsList &lNode = tophits.topHitsLists[iNode];
        for (int64_t iHit = 0; iHit < nCheck && iHit < (int64_t) lNode.hits.size(); iHit++) {
            Besthit bh;
            hitToBestHit(iNode, lNode.hits[iHit], bh);
            setCriterion(nSeqs, bh);
            TopHitsList &lTarget = tophits.topHitsLists[bh.j];

            /* If this criterion is worse than the nCheck-1 entry of the target,
           then skip the check.
           This logic is based on assuming that the list is sorted,
           which is true initially but may not be true later.
           Still, is a good heuristic.
            */
            assert(nCheck > 0);
            assert(nCheck <= (int64_t) lTarget.hits.size());
            Besthit bhCheck;
            hitToBestHit(bh.j, lTarget.hits[nCheck - 1], bhCheck);
            setCriterion(nSeqs, bhCheck);
            if (bhCheck.criterion < bh.criterion) {
                continue;        /* no check needed */
            }

            /* Check if this is present in the top-hit list */
            bool bFound = false;
            for (int64_t iHit2 = 0; iHit2 < (int64_t) lTarget.hits.size() && !bFound; iHit2++) {
                if (lTarget.hits[iHit2].j == iNode) {
                    bFound = true;
                }
            }
            if (!bFound) {
                /* Find the hit with the worst criterion and replace it with this one */
                int64_t iWorst = -1;
                double dWorstCriterion = -1e20;
                for (int64_t iHit2 = 0; iHit2 < (int64_t) lTarget.hits.size(); iHit2++) {
                    Besthit bh2;
                    hitToBestHit(bh.j, lTarget.hits[iHit2], bh2);
                    setCriterion(nSeqs, bh2);
                    if (bh2.criterion > dWorstCriterion) {
                        iWorst = iHit2;
                        dWorstCriterion = bh2.criterion;
                    }
                }
                if (dWorstCriterion > bh.criterion) {
                    assert(iWorst >= 0);
                    lTarget.hits[iWorst].j = iNode;
                    lTarget.hits[iWorst].dist = bh.dist;
                    lReplace++;
                    /* and perhaps update visible */
                    Besthit v;
                    bool bSuccess = getVisible(nSeqs, tophits, bh.j, v);
                    (void) bSuccess;
                    assert(bSuccess);
                    if (bh.criterion < v.criterion) {
                        tophits.visible[bh.j] = lTarget.hits[iWorst];
                    }
                }
            }
        }
    }

    if (options.verbose >= 2) {
        log << strformat("Replaced %ld top hit entries", lReplace) << std::endl;
    }
}


/*
  Find best hit to do in O(N*log(N) + m*L*log(N)) time, by
  copying and sorting the visible list
  updating out-distances for the top (up to m) candidates
  selecting the best hit
  if !fastest then
  	local hill-climbing for a better join,
	using best-hit lists only, and updating
	all out-distances in every best-hit list
*/
AbsNeighbourJoining(void)::topHitNJSearch(int64_t nActive, TopHits &tophits, Besthit &join) {
/* first, do we have at least m/2 candidates in topvisible?
     And remember the best one */
    int64_t nCandidate = 0;
    int64_t iNodeBestCandidate = -1;
    double dBestCriterion = 1e20;

    for (int64_t i = 0; i < (int64_t) tophits.topvisible.size(); i++) {
        int64_t iNode = tophits.topvisible[i];
        Besthit visible;
        if (getVisible(nActive, tophits, iNode, visible)) {
            nCandidate++;
            if (iNodeBestCandidate < 0 || visible.criterion < dBestCriterion) {
                iNodeBestCandidate = iNode;
                dBestCriterion = visible.criterion;
            }
        }
    }

    tophits.topvisibleAge++;
    /* Note we may have only nActive/2 joins b/c we try to store them once */
    if (2 * tophits.topvisibleAge > tophits.m ||
        (3 * nCandidate < (int64_t) tophits.topvisible.size() && 3 * nCandidate < nActive)) {
        /* recompute top visible */
        if (options.verbose > 2) {
            log << strformat("Resetting the top-visible list at nActive=%ld", nActive) << std::endl;
        }

        /* If age is low, then our visible set is becoming too sparse, because we have
           recently recomputed the top visible subset. This is very rare but can happen
           with -fastest. A quick-and-dirty solution is to walk up
           the parents to get additional entries in top hit lists. To ensure that the
           visible set becomes full, pick an arbitrary node if walking up terminates at self.
        */
        if (tophits.topvisibleAge <= 2) {
            if (options.verbose > 2) {
                log << strformat("Expanding visible set by walking up to active nodes at nActive=%ld", nActive)
                    << std::endl;
            }
            for (int64_t iNode = 0; iNode < maxnode; iNode++) {
                if (parent[iNode] >= 0) {
                    continue;
                }
                Hit &v = tophits.visible[iNode];
                int64_t newj = activeAncestor(v.j);
                if (newj >= 0 && newj != v.j) {
                    if (newj == iNode) {
                        /* pick arbitrarily */
                        newj = 0;
                        while (parent[newj] >= 0 || newj == iNode) {
                            newj++;
                        }
                    }
                    assert(newj >= 0 && newj < maxnodes
                           && newj != iNode
                           && parent[newj] < 0);

                    /* Set v to point to newj */
                    Besthit bh = {iNode, newj, (numeric_t) - 1e20, (numeric_t) - 1e20, (numeric_t) - 1e20};
                    setDistCriterion(nActive, bh);
                    v.j = newj;
                    v.dist = bh.dist;
                }
            }
        }
        resetTopVisible(nActive, tophits);
        /* and recurse to try again */
        topHitNJSearch(nActive, tophits, join);
        return;
    }
    if (options.verbose > 2) {
        log << strformat("Top-visible list size %ld (nActive %ld m %ld)", nCandidate, nActive, tophits.m) << std::endl;
    }
    assert(iNodeBestCandidate >= 0 && parent[iNodeBestCandidate] < 0);
    bool bSuccess = getVisible(nActive, tophits, iNodeBestCandidate, join);
    (void) bSuccess;
    assert(bSuccess);
    assert(join.i >= 0 && parent[join.i] < 0);
    assert(join.j >= 0 && parent[join.j] < 0);

    if (options.fastest) {
        return;
    }

    bool changed;
    Besthit join2 = join;
    std::vector<Besthit> bests(options.threads);

    #pragma omp parallel firstprivate(join2) private (changed)
    do {
        changed = false;

        getBestFromTopHits(join2.i, nActive, tophits, bests[omp_get_thread_num()]);
        for (int i = 0; i < options.threads; i++) {
            auto &best = bests[i];
            if (i > 0 && best.i == -1) break;
            assert(best.i == join2.i);
            if (best.j != join2.j && best.criterion < join2.criterion) {
                changed = true;
                if (options.verbose > 2) {
                    log << strformat("BetterI\t%ld\t%ld\t%ld\t%ld\t%f\t%f",
                                     join2.i, join2.j, best.i, best.j,
                                     join2.criterion, best.criterion) << std::endl;
                }
                join2 = best;
            }
        }

        getBestFromTopHits(join2.j, nActive, tophits, bests[omp_get_thread_num()]);
        for (int i = 0; i < options.threads; i++) {
            auto &best = bests[i];
            if (i > 0 && best.i == -1) break;
            assert(best.i == join2.j);
            if (best.j != join2.i && best.criterion < join2.criterion) {
                changed = true;
                if (options.verbose > 2)
                    log << strformat("BetterJ\t%ld\t%ld\t%ld\t%ld\t%f\t%f\n",
                                     join2.i, join2.j, best.i, best.j,
                                     join2.criterion, best.criterion) << std::endl;
                join2 = best;
            }
        }
        if (changed) {
            options.debug.nHillBetter++;
        }
        join = join2;
    } while (changed);
}

/* Updates out-distances but does not reset or update visible set */
AbsNeighbourJoining(void)::getBestFromTopHits(int64_t iNode, int64_t nActive, TopHits &tophits, Besthit &bestjoin) {
    assert(iNode >= 0);
    assert(parent[iNode] < 0);
    TopHitsList &l = tophits.topHitsLists[iNode];
    assert(!l.hits.empty());

    if (!options.fastest) {
        #pragma omp master
        {
            setOutDistance(iNode, nActive); /* ensure out-distances are not stale */
        }
        #pragma omp barrier
    }

    bestjoin.i = -1;
    bestjoin.j = -1;
    bestjoin.dist = (numeric_t) 1e20;
    bestjoin.criterion = (numeric_t) 1e20;

    #pragma omp for schedule(static)
    for (int iBest = 0; iBest < (int64_t) l.hits.size(); iBest++) {
        Besthit bh;
        hitToBestHit(iNode, l.hits[iBest], bh);
        if (updateBestHit(nActive, bh, true)) {
            setCriterion(nActive, bh); /* make sure criterion is correct */
            if (bh.criterion < bestjoin.criterion)
                bestjoin = bh;
        }
    }
    assert(bestjoin.j >= 0);    /* a hit was found */
    assert(bestjoin.i == iNode);
}

/*
  Create a top hit list for the new node, either
  from children (if there are enough best hits left) or by a "refresh"
  Also set visible set for newnode
  Also update visible set for other nodes if we stumble across a "better" hit
*/
AbsNeighbourJoining(void)::topHitJoin(int64_t newnode, int64_t nActive, TopHits &tophits) {//TODO
    int64_t startProfileOps = options.debug.profileOps;
    int64_t startOutProfileOps = options.debug.outprofileOps;
    assert(child[newnode].nChild == 2);
    TopHitsList &lNew = tophits.topHitsLists[newnode];
    assert(lNew.hits.empty());

    /* Copy the hits */
    TopHitsList *lChild[2];
    for (int i = 0; i < 2; i++) {
        lChild[i] = &tophits.topHitsLists[child[newnode].child[i]];
        assert(!lChild[i]->hits.empty());
    }
    int64_t nCombined = lChild[0]->hits.size() + lChild[1]->hits.size();
    std::vector<Besthit> combinedList(nCombined);
    hitsToBestHits(lChild[0]->hits, child[newnode].child[0], combinedList.data());
    hitsToBestHits(lChild[1]->hits, child[newnode].child[1], combinedList.data() + lChild[0]->hits.size());
    /* UniqueBestHits() replaces children (used in the calls to HitsToBestHits)
       with active ancestors, so all distances & criteria will be recomputed */
    std::vector<Besthit> uniqueList;
    uniqueBestHits(nActive, combinedList, uniqueList);
    int64_t nUnique = uniqueList.size();

    combinedList.clear();
    combinedList.reserve(0);

    /* Forget the top-hit lists of the joined nodes */
    for (int i = 0; i < 2; i++) {
        lChild[i]->hits.clear();
    }

    /* Use the average age, rounded up, by 1 Versions 2.0 and earlier
       used the maximum age, which leads to more refreshes without
       improving the accuracy of the NJ phase. Intuitively, if one of
       them was just refreshed then another refresh is unlikely to help.
     */
    lNew.age = (lChild[0]->age + lChild[1]->age + 1) / 2 + 1;

    /* If top hit ages always match (perfectly balanced), then a
       limit of log2(m) would mean a refresh after
       m joins, which is about what we want.
    */
    int64_t tophitAgeLimit = std::max((int64_t) 1, (int64_t) (0.5 + std::log((double) tophits.m) / std::log(2.0)));

    /* Either use the merged list as candidate top hits, or
       move from 2nd level to 1st level, or do a refresh
       UniqueBestHits eliminates hits to self, so if nUnique==nActive-1,
       we've already done the exhaustive search.

       Either way, we set tophits, visible(newnode), update visible of its top hits,
       and modify topvisible: if we do a refresh, then we reset it, otherwise we update
    */
    bool bSecondLevel = lChild[0]->hitSource >= 0 && lChild[1]->hitSource >= 0;
    bool bUseUnique = nUnique == nActive - 1
                      || (lNew.age <= tophitAgeLimit
                          && nUnique >= (bSecondLevel ? (int64_t) (0.5 + options.tophits2Refresh * tophits.q)
                                                      : (int64_t) (0.5 + tophits.m * options.tophitsRefresh)));
    if (bUseUnique && options.verbose > 2) {
        log << strformat("Top hits for %ld from combined %ld nActive=%ld tophitsage %ld %s",
                         newnode, nUnique, nActive, lNew.age, bSecondLevel ? "2ndlevel" : "1stlevel") << std::endl;
    }

    if (!bUseUnique && bSecondLevel && lNew.age <= tophitAgeLimit) {
        int64_t source = activeAncestor(lChild[0]->hitSource);
        if (source == newnode) {
            source = activeAncestor(lChild[1]->hitSource);
        }
        /* In parallel mode, it is possible that we would select a node as the
           hit-source and then over-write that top hit with a short list.
           So we need this sanity check.
        */
        if (source != newnode && source >= 0 && tophits.topHitsLists[source].hitSource < 0) {

            /* switch from 2nd-level to 1st-level top hits -- compute top hits list
           of node from what we have so far plus the active source plus its top hits */
            TopHitsList &lSource = tophits.topHitsLists[source];
            assert(lSource.hitSource < 0);
            assert(!lSource.hits.empty());
            int64_t nMerge = (int64_t) lSource.hits.size() + nUnique + 1;
            std::vector<Besthit> mergeList(uniqueList); /* Copy */
            mergeList.resize(nMerge);

            int64_t iMerge = nUnique;
            mergeList[iMerge].i = newnode;
            mergeList[iMerge].j = source;
            setDistCriterion(nActive, mergeList[iMerge]);
            iMerge++;
            hitsToBestHits(lSource.hits, newnode, mergeList.data() + iMerge);
            for (int64_t i = 0; i < (int64_t) lSource.hits.size(); i++) {
                setDistCriterion(nActive, mergeList[iMerge]);
                iMerge++;
            }
            assert(iMerge == nMerge);

            uniqueList.clear();

            uniqueBestHits(nActive, mergeList, uniqueList);

            mergeList.clear();
            mergeList.reserve(0);

            assert(nUnique > 0);
            bUseUnique = nUnique >= (int64_t) (0.5 + tophits.m * options.tophitsRefresh);
            bSecondLevel = false;

            if (bUseUnique && options.verbose > 2) {
                log << strformat("Top hits for %ld from children and source %ld's %zd hits, nUnique %ld",
                                 newnode, source, lSource.hits.size(), nUnique);
            }
        }
    }

    if (bUseUnique) {
        if (bSecondLevel) {
            /* pick arbitrarily */
            lNew.hitSource = lChild[0]->hitSource;
        }
        int64_t nSave = std::min(nUnique, bSecondLevel ? tophits.q : tophits.m);
        assert(nSave > 0);
        if (options.verbose > 2 && options.threads == 1) {
            log << strformat("Combined %ld ops so far %ld\n", nUnique, options.debug.profileOps - startProfileOps)
                << std::endl;
        }
        sortSaveBestHits(newnode, uniqueList, nUnique, nSave, tophits);
        assert(!lNew.hits.empty()); /* set by sort/save */
        tophits.visible[newnode] = lNew.hits[0];
        updateTopVisible(nActive, newnode, tophits.visible[newnode], tophits);
        uniqueList.resize(nSave);
        updateVisible(nActive, uniqueList, tophits);
    } else {
        /* need to refresh: set top hits for node and for its top hits */
        if (options.verbose > 2) {
            log << strformat("Top hits for %ld by refresh (%ld unique age %ld) nActive=%ld",
                             newnode, nUnique, lNew.age, nActive) << std::endl;
        }
        options.debug.nRefreshTopHits++;
        lNew.age = 0;

        /* ensure all out-distances are up to date ahead of time
           to avoid any data overwriting issues.
        */
        #pragma omp parallel for schedule(dynamic)
        for (int64_t iNode = 0; iNode < maxnode; iNode++) {
            if (parent[iNode] < 0) {
                if (options.fastest) {
                    Besthit bh;
                    bh.i = iNode;
                    bh.j = iNode;
                    bh.dist = 0;
                    setCriterion(nActive, bh);
                } else {
                    setOutDistance(iNode, nActive);
                }
            }
        }

        /* exhaustively get the best 2*m hits for newnode, set visible, and save the top m */
        std::vector<Besthit> allhits(maxnode);
        assert(2 * tophits.m <= maxnode);
        Besthit bh;
        setBestHit(newnode, nActive, bh, allhits.data());
        psort(allhits.begin(), allhits.end(), CompareHitsByCriterion());
        sortSaveBestHits(newnode, allhits, maxnode, tophits.m, tophits);

        /* Do not need to call UpdateVisible because we set visible below */

        /* And use the top 2*m entries to expand other best-hit lists, but only for top m */
        #pragma omp parallel for schedule(dynamic)
        for (int64_t iHit = 0; iHit < tophits.m; iHit++) {
            if (allhits[iHit].i < 0) {
                continue;
            }
            int64_t iNode = allhits[iHit].j;
            assert(iNode >= 0);
            if (parent[iNode] >= 0) {
                continue;
            }
            TopHitsList &l = tophits.topHitsLists[iNode];
            int64_t nHitsOld = l.hits.size();
            assert(nHitsOld <= tophits.m);
            l.age = 0;

            /* Merge: old hits into 0.nHitsOld and hits from iNode above that */
            std::vector<Besthit> bothList(3 * tophits.m);
            hitsToBestHits(l.hits, iNode, bothList.data()); /* does not compute criterion */
            for (int64_t i = 0; i < nHitsOld; i++) {
                setCriterion(nActive, bothList[i]);
            }
            if (nActive <= 2 * tophits.m) {
                l.hitSource = -1;    /* abandon the 2nd-level top-hits heuristic */
            }
            int64_t nNewHits = l.hitSource >= 0 ? tophits.q : tophits.m;
            assert(nNewHits > 0);

            transferBestHits(nActive, iNode, allhits, 2 * nNewHits, &bothList[nHitsOld], false);
            /* rely on UniqueBestHits to update dist and/or criterion */
            std::vector<Besthit> uniqueList2;
            bothList.resize(nHitsOld + 2 * nNewHits);
            uniqueBestHits(nActive, bothList, uniqueList2);
            assert(!uniqueList2.empty());

            /* Note this will overwrite l, but we saved nHitsOld */
            sortSaveBestHits(iNode, uniqueList2, uniqueList2.size(), nNewHits, tophits);
            /* will update topvisible below */
            tophits.visible[iNode] = tophits.topHitsLists[iNode].hits[0];
        }

        resetTopVisible(nActive, tophits); /* outside of the parallel phase */
    }
    if (options.verbose > 2) {
        log << "New top-hit list for " << newnode;
        if (options.threads == 1) {
            log << strformat("profile-ops %ld (out-ops %ld)",
                             options.debug.profileOps - startProfileOps,
                             options.debug.outprofileOps - startOutProfileOps);
        }
        log << strformat(": source %ld age %ld members ",
                         lNew.hitSource, lNew.age);
        for (int64_t i = 0; i < (int64_t) lNew.hits.size(); i++) {
            log << " " << lNew.hits[i].j;
        }
        log << std::endl;
    }
}

AbsNeighbourJoining(void)::sortSaveBestHits(int64_t iNode, std::vector<Besthit> &besthits, int64_t nIn, int64_t nOut,
                                            TopHits &tophits, bool sort) {
    assert(nIn > 0);
    assert(nOut > 0);

    if (sort) {
        psort(besthits.begin(), besthits.end(), CompareHitsByCriterion());
    }

    /* First count how many we will save
       Not sure if removing duplicates is actually necessary.
     */
    int64_t nSave = 0;
    int64_t jLast = -1;
    for (int64_t iBest = 0; iBest < nIn && nSave < nOut; iBest++) {
        if (besthits[iBest].i < 0) {
            continue;
        }
        assert(besthits[iBest].i == iNode);
        int64_t j = besthits[iBest].j;
        if (j != iNode && j != jLast && j >= 0) {
            nSave++;
            jLast = j;
        }
    }

    assert(nSave > 0);

    TopHitsList &l = tophits.topHitsLists[iNode];
    l.hits.resize(nSave);

    int64_t iSave = 0;
    jLast = -1;
    for (int64_t iBest = 0; iBest < nIn && iSave < nSave; iBest++) {
        int64_t j = besthits[iBest].j;
        if (j != iNode && j != jLast && j >= 0) {
            l.hits[iSave].j = j;
            l.hits[iSave].dist = besthits[iBest].dist;
            iSave++;
            jLast = j;
        }
    }
    assert(iSave == nSave);
}

AbsNeighbourJoining(void)::transferBestHits(int64_t nActive, int64_t iNode, std::vector<Besthit> &oldhits,
                                            int64_t nOldHits, Besthit newhits[], bool updateDistances) {
    assert(iNode >= 0);
    assert(parent[iNode] < 0);

    for (int64_t iBest = 0; iBest < nOldHits; iBest++) {
        Besthit &oldhit = oldhits[iBest];
        Besthit &newhit = newhits[iBest];
        newhit.i = iNode;
        newhit.j = activeAncestor(oldhit.j);
        newhit.dist = oldhit.dist;    /* may get reset below */
        newhit.weight = oldhit.weight;
        newhit.criterion = oldhit.criterion;

        if (newhit.j < 0 || newhit.j == iNode) {
            newhit.weight = 0;
            newhit.dist = (numeric_t) - 1e20;
            newhit.criterion = (numeric_t) 1e20;
        } else if (newhit.i != oldhit.i || newhit.j != oldhit.j) {
            if (updateDistances) {
                setDistCriterion(nActive, newhit);
            } else {
                newhit.dist = (numeric_t) - 1e20;
                newhit.criterion = (numeric_t) 1e20;
            }
        } else {
            if (updateDistances) {
                setCriterion(nActive, newhit);
            } else {
                newhit.criterion = (numeric_t) 1e20;    /* leave dist alone */
            }
        }
    }
}

AbsNeighbourJoining(void)::hitsToBestHits(std::vector<Hit> &hits, int64_t iNode, Besthit newhits[]) {
    for (int64_t i = 0; i < (int64_t) hits.size(); i++) {
        Hit &hit = hits[i];
        Besthit &bh = newhits[i];
        bh.i = iNode;
        bh.j = hit.j;
        bh.dist = hit.dist;
        bh.criterion = (numeric_t) 1e20;
        bh.weight = -1;        /* not the true value -- we compute these directly when needed */
    }
}

AbsNeighbourJoining(void)::hitToBestHit(int64_t i, Hit &hit, Besthit &out) {
    out.i = i;
    out.j = hit.j;
    out.dist = hit.dist;
    out.criterion = (numeric_t) 1e20;
    out.weight = -1;
}

AbsNeighbourJoining(void)::updateVisible(int64_t nActive, std::vector<Besthit> &tophitsNode, TopHits &tophits) {
    for (int64_t iHit = 0; iHit < (int64_t) tophitsNode.size(); iHit++) {
        Besthit &hit = tophitsNode[iHit];
        if (hit.i < 0) {
            continue;
        }    /* possible empty entries */
        assert(parent[hit.i] < 0);
        assert(hit.j >= 0 && parent[hit.j] < 0);
        Besthit visible;
        bool bSuccess = getVisible(nActive, tophits, hit.j, visible);
        if (!bSuccess || hit.criterion < visible.criterion) {
            if (bSuccess) {
                options.debug.nVisibleUpdate++;
            }
            Hit &v = tophits.visible[hit.j];
            v.j = hit.i;
            v.dist = hit.dist;
            updateTopVisible(nActive, hit.j, v, tophits);
            if (options.verbose > 5) {
                log << strformat("NewVisible %ld %ld %f", hit.j, v.j, v.dist) << std::endl;
            }
        }
    } /* end loop over hits */
}

/* Update the top-visible list to perhaps include visible[iNode] */
AbsNeighbourJoining(void)::updateTopVisible(int64_t nActive, int64_t iIn, Hit &hit, TopHits &tophits) {
    bool bIn = false;        /* placed in the list */

    /* First, if the list is not full, put it in somewhere */
    for (int64_t i = 0; i < (int64_t) tophits.topvisible.size() && !bIn; i++) {
        int64_t iNode = tophits.topvisible[i];
        if (iNode == iIn) {
            /* this node is already in the top hit list */
            bIn = true;
        } else if (iNode < 0 || parent[iNode] >= 0) {
            /* found an empty spot */
            bIn = true;
            tophits.topvisible[i] = iIn;
        }
    }

    int64_t iPosWorst = -1;
    double dCriterionWorst = -1e20;
    if (!bIn) {
        /* Search for the worst hit */
        for (int64_t i = 0; i < (int64_t) tophits.topvisible.size() && !bIn; i++) {
            int64_t iNode = tophits.topvisible[i];
            assert(iNode >= 0 && parent[iNode] < 0 && iNode != iIn);
            Besthit visible;
            if (!getVisible(nActive, tophits, iNode, visible)) {
                /* found an empty spot */
                tophits.topvisible[i] = iIn;
                bIn = true;
            } else if (visible.i == hit.j && visible.j == iIn) {
                /* the reverse hit is already in the top hit list */
                bIn = true;
            } else if (visible.criterion >= dCriterionWorst) {
                iPosWorst = i;
                dCriterionWorst = visible.criterion;
            }
        }
    }

    if (!bIn && iPosWorst >= 0) {
        Besthit visible;
        hitToBestHit(iIn, hit, visible);
        setCriterion(nActive, visible);
        if (visible.criterion < dCriterionWorst) {
            if (options.verbose > 2) {
                int64_t iOld = tophits.topvisible[iPosWorst];
                log << strformat("TopVisible replace %ld=>%ld with %ld=>%ld",
                                 iOld, tophits.visible[iOld].j, visible.i, visible.j) << std::endl;
            }
            tophits.topvisible[iPosWorst] = iIn;
        }
    }

    if (options.verbose > 2) {
        log << "Updated TopVisible: ";
        for (int64_t i = 0; i < (int64_t) tophits.topvisible.size(); i++) {
            int iNode = tophits.topvisible[i];
            if (iNode >= 0 && parent[iNode] < 0) {
                Besthit bh;
                hitToBestHit(iNode, tophits.visible[iNode], bh);
                setDistCriterion(nActive, bh);
                log << strformat(" %ld=>%ld:%.4f", bh.i, bh.j, bh.criterion);
            }
        }
        log << std::endl;
    }
}

AbsNeighbourJoining(void)::resetTopVisible(int64_t nActive, TopHits &tophits) {
    std::vector<Besthit> visibleSorted(nActive);
    int64_t nVisible = 0;        /* #entries in visibleSorted */
    for (int64_t iNode = 0; iNode < maxnode; iNode++) {
        /* skip joins involving stale nodes */
        if (parent[iNode] >= 0) {
            continue;
        }
        Besthit v;
        if (getVisible(nActive, tophits, iNode, v)) {
            assert(nVisible < nActive);
            visibleSorted[nVisible++] = v;
        }
    }
    assert(nVisible > 0);

    psort(visibleSorted.begin(), visibleSorted.end(), CompareHitsByCriterion());

    /* Only keep the top m items, and try to avoid duplicating i->j with j->i
       Note that visible(i) -> j does not necessarily imply visible(j) -> i,
       so we store what the pairing was (or -1 for not used yet)
     */
    std::vector<int64_t> inTopVisible(maxnodes);
    for (int64_t i = 0; i < maxnodes; i++)
        inTopVisible[i] = -1;

    if (options.verbose > 2) {
        log << strformat("top-hit search: nActive %ld nVisible %ld considering up to %ld items",
                         nActive, nVisible, tophits.m) << std::endl;
    }

    /* save the sorted indices in topvisible */
    int64_t iSave = 0;
    for (int64_t i = 0; i < nVisible && iSave < (int64_t) tophits.topvisible.size(); i++) {
        Besthit &v = visibleSorted[i];
        if (inTopVisible[v.i] != v.j) { /* not seen already */
            tophits.topvisible[iSave++] = v.i;
            inTopVisible[v.i] = v.j;
            inTopVisible[v.j] = v.i;
        }
    }
    while (iSave < (int64_t) tophits.topvisible.size()) {
        tophits.topvisible[iSave++] = -1;
    }
    tophits.topvisibleAge = 0;
    if (options.verbose > 2) {
        log << "Reset TopVisible: ";
        for (int64_t i = 0; i < (int64_t) tophits.topvisible.size(); i++) {
            int64_t iNode = tophits.topvisible[i];
            if (iNode < 0) {
                break;
            }
            log << strformat(" %ld=>%ld", iNode, tophits.visible[iNode].j);
        }
        log << std::endl;;
    }
}

AbsNeighbourJoining(void)::uniqueBestHits(int64_t nActive, std::vector<Besthit> &combined, std::vector<Besthit> &out) {
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (int64_t iHit = 0; iHit < (int64_t) combined.size(); iHit++) {
            Besthit &hit = combined[iHit];
            updateBestHit(nActive, hit, false);
        }

        #pragma omp master
        {
            psort(combined.begin(), combined.end(), CompareHitsByIJ(), true);

            out.reserve(combined.size());
            int64_t iSavedLast = -1;

            /* First build the new list */
            for (int64_t iHit = 0; iHit < (int64_t) combined.size(); iHit++) {
                Besthit &hit = combined[iHit];
                if (hit.i < 0 || hit.j < 0) {
                    continue;
                }
                if (iSavedLast >= 0) {
                    /* toss out duplicates */
                    Besthit &saved = combined[iSavedLast];
                    if (saved.i == hit.i && saved.j == hit.j) {
                        continue;
                    }
                }
                assert(hit.j >= 0 && parent[hit.j] < 0);
                out.push_back(hit);
                iSavedLast = iHit;
            }
        }
        #pragma omp barrier

        /* Then do any updates to the criterion or the distances in parallel */
        #pragma omp for schedule(dynamic)
        for (int64_t iHit = 0; iHit < (int64_t) out.size(); iHit++) {
            Besthit &hit = out[iHit];
            if (hit.dist < 0.0) {
                setDistCriterion(nActive, hit);
            } else {
                setCriterion(nActive, hit);
            }
        }
    }
}


AbsNeighbourJoining(typename veryfasttree::NeighbourJoining<Precision, Operations>::NNI)::chooseNNI(
        Profile *profiles4[4],
        double criteria[3]) {
    double d[6];
    correctedPairDistances(profiles4, 4, d);
    double penalty[3];        /* indexed as nni_t */
    quartetConstraintPenalties(profiles4, penalty);
    criteria[ABvsCD] = d[qAB] + d[qCD] + penalty[ABvsCD];
    criteria[ACvsBD] = d[qAC] + d[qBD] + penalty[ACvsBD];
    criteria[ADvsBC] = d[qAD] + d[qBC] + penalty[ADvsBC];

    auto choice = ABvsCD;
    if (criteria[ACvsBD] < criteria[ABvsCD] && criteria[ACvsBD] <= criteria[ADvsBC]) {
        choice = ACvsBD;
    } else if (criteria[ADvsBC] < criteria[ABvsCD] && criteria[ADvsBC] <= criteria[ACvsBD]) {
        choice = ADvsBC;
    }
    if (options.verbose > 1 && penalty[choice] > penalty[ABvsCD] + 1e-6) {
        log << strformat("Worsen constraint: from %.3f to %.3f distance %.3f to %.3f: ",
                         penalty[ABvsCD], penalty[choice],
                         criteria[ABvsCD], choice == ACvsBD ? criteria[ACvsBD] : criteria[ADvsBC]);

        for (int64_t iC = 0; iC < (int64_t) constraintSeqs.size(); iC++) {
            double ppart[3];
            if (quartetConstraintPenaltiesPiece(profiles4, iC, ppart)) {
                double old_penalty = ppart[ABvsCD];
                double new_penalty = ppart[choice];
                if (new_penalty > old_penalty + 1e-6) {
                    log << strformat(" %ld (%ld/%ld %ld/%ld %ld/%ld %ld/%ld)", iC,
                                     profiles4[0]->nOn[iC], profiles4[0]->nOff[iC],
                                     profiles4[1]->nOn[iC], profiles4[1]->nOff[iC],
                                     profiles4[2]->nOn[iC], profiles4[2]->nOff[iC],
                                     profiles4[3]->nOn[iC], profiles4[3]->nOff[iC]);
                }
            }
        }
        log << std::endl;
    }
    if (options.verbose > 3) {
        log << strformat("NNI scores ABvsCD %.5f ACvsBD %.5f ADvsBC %.5f choice %s",
                         criteria[ABvsCD], criteria[ACvsBD], criteria[ADvsBC],
                         choice == ABvsCD ? "AB|CD" : (choice == ACvsBD ? "AC|BD" : "AD|BC")) << std::endl;
    }
    return (choice);
}

AbsNeighbourJoining(typename veryfasttree::NeighbourJoining<Precision, Operations>::NNI)::
MLQuartetNNI(Profile *profiles4[4], double criteria[3], numeric_t len[5], bool bFast) {
    double lenABvsCD[5] = {len[LEN_A], len[LEN_B], len[LEN_C], len[LEN_D], len[LEN_I]};
    double lenACvsBD[5] = {len[LEN_A], len[LEN_C], len[LEN_B], len[LEN_D], len[LEN_I]};   /* Swap B & C */
    double lenADvsBC[5] = {len[LEN_A], len[LEN_D], len[LEN_C], len[LEN_B], len[LEN_I]};   /* Swap B & D */
    bool bConsiderAC = true;
    bool bConsiderAD = true;

    int64_t nRounds = options.mlAccuracy < 2 ? 2 : options.mlAccuracy;
    double penalty[3];
    quartetConstraintPenalties(profiles4, penalty);
    if (penalty[ABvsCD] > penalty[ACvsBD] || penalty[ABvsCD] > penalty[ADvsBC]) {
        bFast = false; /* turn off star topology test */
    }

    int64_t iRound;
    for (iRound = 0; iRound < nRounds; iRound++) {
        if (omp_in_parallel() || options.threads == 1) {
            bool bStarTest = false;
            criteria[ABvsCD] = MLQuartetOptimize(*profiles4[0], *profiles4[1], *profiles4[2], *profiles4[3],
                                                 lenABvsCD, &bStarTest,  /*site_likelihoods*/nullptr)
                               - penalty[ABvsCD];    /* subtract penalty b/c we are trying to maximize log lk */
            if (bStarTest) {
                options.debug.nStarTests++;
                criteria[ACvsBD] = -1e20;
                criteria[ADvsBC] = -1e20;
                len[LEN_I] = lenABvsCD[LEN_I];
                return ABvsCD;
            }

            if (bConsiderAC) {
                criteria[ACvsBD] = MLQuartetOptimize(*profiles4[0], *profiles4[2], *profiles4[1], *profiles4[3],
                                                     lenACvsBD, nullptr, /*site_likelihoods*/nullptr)
                                   - penalty[ACvsBD];
            }

            if (bConsiderAD) {
                criteria[ADvsBC] = MLQuartetOptimize(*profiles4[0], *profiles4[3], *profiles4[2], *profiles4[1],
                                                     lenADvsBC, nullptr, /*site_likelihoods*/nullptr)
                                   - penalty[ADvsBC];
            }
        } else {
            #pragma omp parallel
            #pragma omp sections
            {
                #pragma omp section
                {
                    criteria[ABvsCD] = MLQuartetOptimize(*profiles4[0], *profiles4[1], *profiles4[2], *profiles4[3],
                                                         lenABvsCD, nullptr,  /*site_likelihoods*/nullptr)
                                       - penalty[ABvsCD];    /* subtract penalty b/c we are trying to maximize log lk */
                }

                #pragma omp section
                {
                    if (bConsiderAC) {
                        criteria[ACvsBD] = MLQuartetOptimize(*profiles4[0], *profiles4[2], *profiles4[1], *profiles4[3],
                                                             lenACvsBD, nullptr, /*site_likelihoods*/nullptr)
                                           - penalty[ACvsBD];
                    }
                }

                #pragma omp section
                {
                    if (bConsiderAD) {
                        criteria[ADvsBC] = MLQuartetOptimize(*profiles4[0], *profiles4[3], *profiles4[2], *profiles4[1],
                                                             lenADvsBC, nullptr, /*site_likelihoods*/nullptr)
                                           - penalty[ADvsBC];
                    }
                }
            }
        }/* end parallel sections */

        if (options.mlAccuracy < 2) {
            /* If clearly worse then ABvsCD, or have short internal branch length and worse, then
               give up */
            if (criteria[ACvsBD] < criteria[ABvsCD] - Constants::closeLogLkLimit
                || (lenACvsBD[LEN_I] <= 2.0 * options.MLMinBranchLength && criteria[ACvsBD] < criteria[ABvsCD])) {
                bConsiderAC = false;
            }
            if (criteria[ADvsBC] < criteria[ABvsCD] - Constants::closeLogLkLimit
                || (lenADvsBC[LEN_I] <= 2.0 * options.MLMinBranchLength && criteria[ADvsBC] < criteria[ABvsCD])) {
                bConsiderAD = false;
            }
            if (!bConsiderAC && !bConsiderAD) {
                break;
            }
            /* If clearly better than either alternative, then give up
               (Comparison is probably biased in favor of ABvsCD anyway) */
            if (criteria[ACvsBD] > criteria[ABvsCD] + Constants::closeLogLkLimit
                && criteria[ACvsBD] > criteria[ADvsBC] + Constants::closeLogLkLimit) {
                break;
            }
            if (criteria[ADvsBC] > criteria[ABvsCD] + Constants::closeLogLkLimit
                && criteria[ADvsBC] > criteria[ACvsBD] + Constants::closeLogLkLimit) {
                break;
            }
        }
    } /* end loop over rounds */

    if (options.verbose > 2) {
        log << strformat("Optimized quartet for %ld rounds: ABvsCD %.5f ACvsBD %.5f ADvsBC %.5f",
                         iRound, criteria[ABvsCD], criteria[ACvsBD], criteria[ADvsBC]) << std::endl;
    }
    if (criteria[ACvsBD] > criteria[ABvsCD] && criteria[ACvsBD] > criteria[ADvsBC]) {
        for (int i = 0; i < 5; i++) {
            len[i] = lenACvsBD[i];
        }
        return (ACvsBD);
    } else if (criteria[ADvsBC] > criteria[ABvsCD] && criteria[ADvsBC] > criteria[ACvsBD]) {
        for (int i = 0; i < 5; i++) {
            len[i] = lenADvsBC[i];
        }
        return (ADvsBC);
    } else {
        for (int i = 0; i < 5; i++) {
            len[i] = lenABvsCD[i];
        }
        return (ABvsCD);
    }
}

AbsNeighbourJoining(inline int64_t)::traverseOptimizeAllBranchLengths(int64_t node,
                                                                      std::unique_ptr<Profile> upProfiles[],
                                                                      std::vector<ibool> &traversal) {
    const bool parallel = omp_in_parallel();
    int64_t branchRoot = node;
    int64_t iNodesDone = 0;
    while ((node = traversePostorder(node,/*IN/OUT*/traversal, /*pUp*/nullptr, branchRoot)) >= 0) {//level-3
        int64_t nChild = child[node].nChild;
        if (nChild > 0) {
            if (iNodesDone > 0 && (iNodesDone % 100) == 0 && !parallel) {
                progressReport.print("ML Lengths %ld of %ld splits", iNodesDone + 1, (int64_t) (maxnode - nSeqs));
            }
            iNodesDone++;

            /* optimize the branch lengths between self, parent, and children,
               with two iterations
            */
            assert(nChild == 2 || nChild == 3);
            int64_t nodes[3] = {child[node].child[0],
                                child[node].child[1],
                                nChild == 3 ? child[node].child[2] : node};
            Profile *profiles3[3] = {&profiles[nodes[0]],
                                     &profiles[nodes[1]],
                                     nChild == 3 ? &profiles[nodes[2]]
                                                 : getUpProfile(/*IN/OUT*/upProfiles, node, /*useML*/true)};
            for (int iter = 0; iter < 2; iter++) {
                for (int i = 0; i < 3; i++) {
                    Profile &pA = *profiles3[i];
                    int b1 = (i + 1) % 3;
                    int b2 = (i + 2) % 3;
                    Profile pB(nPos, constraintSeqs.size());
                    posteriorProfile(pB, *profiles3[b1], *profiles3[b2], branchlength[nodes[b1]],
                                     branchlength[nodes[b2]]);
                    double len = branchlength[nodes[i]];
                    if (len < options.MLMinBranchLength) {
                        len = options.MLMinBranchLength;
                    }
                    MLPairOptimize(pA, pB, /*IN/OUT*/&len);
                    branchlength[nodes[i]] = len;
                    if (options.verbose > 3) {
                        log << strformat("Optimize length for %ld to %.3f", nodes[i], branchlength[nodes[i]])
                            << std::endl;
                    }
                }
            }
            if (node != root) {
                recomputeProfile(/*IN/OUT*/upProfiles, node, /*useML*/true);
                upProfiles[node].reset();
            }
        }
    }
    return iNodesDone;
}

AbsNeighbourJoining(void)::optimizeAllBranchLengths() {
    if (nSeqs < 2) {
        return;
    }
    if (nSeqs == 2) {
        int64_t parent = root;
        assert(child[parent].nChild == 2);
        int64_t nodes[2] = {child[parent].child[0], child[parent].child[1]};
        double length = 1.0;
        MLPairOptimize(profiles[nodes[0]], profiles[nodes[1]], /*IN/OUT*/&length);
        branchlength[nodes[0]] = length / 2.0;
        branchlength[nodes[1]] = length / 2.0;
        return;
    };

    std::vector<ibool> traversal(maxnodes, false);

    if (options.threads > 1 && options.threadsLevel > 2) {
        std::vector<std::vector<int64_t>> chunks;
        treePartition(chunks, traversal, 1);
        int64_t done = 0;

        #pragma omp parallel
        {
            std::vector<std::unique_ptr<Profile>> upProfiles2(maxnodes);
            int64_t done2;

            #pragma omp for schedule(static, 1)
            for (int i = 0; i < (int) chunks.size(); i++) {
                for (int64_t subtreeRoot: chunks[i]) {
                    for (int i = 0; i < child[subtreeRoot].nChild; i++) {
                        done2 = traverseOptimizeAllBranchLengths(child[subtreeRoot].child[i],
                                                                 upProfiles2.data(), traversal);
                        if (options.verbose > 0) {
                            #pragma omp critical
                            {
                                done += done2;
                                progressReport.print("ML Lengths %ld of %ld splits", done + 1,
                                                     (int64_t) (maxnode - nSeqs));
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<std::unique_ptr<Profile>> upProfiles(maxnodes);
    traverseOptimizeAllBranchLengths(root, upProfiles.data(), traversal);
}

AbsNeighbourJoining(inline double)::traverseTreeLogLk(int64_t node, std::vector<double> &site_likelihood,
                                                      double site_loglk[], std::vector<ibool> &traversal) {
    int64_t branchRoot = node;
    double loglk = 0;
    while ((node = traversePostorder(node, traversal, /*pUp*/nullptr, branchRoot)) >= 0) {//level-1
        int64_t nChild = child[node].nChild;
        if (nChild == 0) {
            continue;
        }
        assert(nChild >= 2);
        int64_t *children = child[node].child;
        double loglkchild = pairLogLk(profiles[children[0]], profiles[children[1]],
                                      branchlength[children[0]] + branchlength[children[1]], site_likelihood.data());
        loglk += loglkchild;
        if (!site_likelihood.empty()) {
            /* prevent underflows */
            for (int64_t i = 0; i < nPos; i++) {
                while (site_likelihood[i] < Constants::LkUnderflow) {
                    site_likelihood[i] *= Constants::LkUnderflowInv;
                    site_loglk[i] -= Constants::LogLkUnderflow;
                }
            }
        }
        if (options.verbose > 2) {
            log << strformat("At %ld: LogLk(%ld:%.4f,%ld:%.4f) = %.3f",
                             node,
                             children[0], branchlength[children[0]],
                             children[1], branchlength[children[1]],
                             loglkchild) << std::endl;
        }
        if (child[node].nChild == 3) {
            assert(node == root);
            /* Infer the common parent of the 1st two to define the third... */
            Profile AB(nPos,/*nConstraints*/0);
            posteriorProfile(AB, profiles[children[0]], profiles[children[1]], branchlength[children[0]],
                             branchlength[children[1]]);
            double loglkup = pairLogLk(AB, profiles[children[2]], branchlength[children[2]], site_likelihood.data());
            loglk += loglkup;
            if (options.verbose > 2) {
                log << strformat("At root %ld: LogLk((%ld/%ld),%ld:%.3f) = %.3f",
                                 node, children[0], children[1], children[2],
                                 branchlength[children[2]], loglkup) << std::endl;
            }
        }
    }
    return loglk;
}

AbsNeighbourJoining(double)::treeLogLk(double site_loglk[]) {
    if (nSeqs < 2) {
        return 0.0;
    }
    double loglk = 0.0;
    std::vector<double> site_likelihood;

    if (site_loglk != nullptr) {
        site_likelihood.resize(nPos, 1.0);

        for (int64_t i = 0; i < nPos; i++) {
            site_likelihood[i] = 1.0;
            site_loglk[i] = 0.0;
        }
    }

    std::vector<ibool> traversal(maxnodes, false);

    if (options.threads > 1 && options.threadsLevel > 2) {
        std::vector<std::vector<int64_t>> chunks;
        treePartition(chunks, traversal, 0);

        #pragma omp parallel
        {
            double loglk2 = 0.0;
            std::vector<double> site_likelihood2 = site_likelihood;
            std::vector<double> site_loglk2(site_loglk != nullptr ? nPos : 0, 0.0);

            #pragma omp for schedule(static, 1)
            for (int i = 0; i < (int) chunks.size(); i++) {
                for (int64_t subtreeRoot: chunks[i]) {
                    loglk2 += traverseTreeLogLk(subtreeRoot, site_likelihood2, site_loglk2.data(), traversal);
                }
            }
            #pragma omp critical
            {
                loglk += loglk2;
                if (!site_likelihood.empty()) {
                    for (int64_t i = 0; i < nPos; i++) {
                        site_likelihood[i] *= site_likelihood2[i];
                        site_loglk[i] += site_loglk2[i];
                    }
                }
            }
        }

    }

    loglk += traverseTreeLogLk(root, site_likelihood, site_loglk, traversal);

    if (!site_likelihood.empty()) {
        for (int64_t i = 0; i < nPos; i++) {
            site_loglk[i] += std::log(site_likelihood[i]);
        }
    }

    /* For Jukes-Cantor, with a tree of size 4, if the children of the root are
       (A,B), C, and D, then
       P(ABCD) = P(A) P(B|A) P(C|AB) P(D|ABC)

       Above we compute P(B|A) P(C|AB) P(D|ABC) -- note P(B|A) is at the child of root
       and P(C|AB) P(D|ABC) is at root.

       Similarly if the children of the root are C, D, and (A,B), then
       P(ABCD) = P(C|D) P(A|B) P(AB|CD) P(D), and above we compute that except for P(D)

       So we need to multiply by P(A) = 0.25, so we pay log(4) at each position
       (if ungapped). Each gapped position in any sequence reduces the payment by log(4)

       For JTT or GTR, we are computing P(A & B) and the posterior profiles are scaled to take
       the prior into account, so we do not need any correction.
       codeFreq[NOCODE] is scaled x higher so that P(-) = 1 not P(-)=1/nCodes, so gaps
       do not need to be corrected either.
     */

    if (options.nCodes == 4 && !transmat) {
        int64_t nGaps = 0;
        double logNCodes = std::log((double) options.nCodes);
        for (int64_t i = 0; i < nPos; i++) {
            int nGapsThisPos = 0;
            for (int64_t node = 0; node < (int64_t) nSeqs; node++) {
                auto &codes = profiles[node].codes;
                if (codes[i] == NOCODE) {
                    nGapsThisPos++;
                }
            }
            nGaps += nGapsThisPos;
            if (site_loglk != nullptr) {
                site_loglk[i] += nGapsThisPos * logNCodes;
                if (options.nCodes == 4 && !transmat) {
                    site_loglk[i] -= logNCodes;
                }
            }
        }
        loglk -= nPos * logNCodes;
        loglk += nGaps * logNCodes;    /* do not pay for gaps -- only Jukes-Cantor */
    }
    return loglk;
}

AbsNeighbourJoining(double)::gammaLogLk(Siteratelk &s, double gamma_loglk_sites[]) {
    std::vector<double> dRate(options.nRateCats);
    for (int64_t iRate = 0; iRate < (int64_t) dRate.size(); iRate++) {
        /* The probability density for each rate is approximated by the total
           density between the midpoints */
        double pMin = iRate == 0 ? 0.0 : pGamma(s.mult * (s.rates[iRate - 1] + s.rates[iRate]) / 2.0, s.alpha);
        double pMax = iRate == (int64_t) (rates.rates.size() - 1) ? 1.0 :
                      pGamma(s.mult * (s.rates[iRate] + s.rates[iRate + 1]) / 2.0, s.alpha);
        dRate[iRate] = pMax - pMin;
    }

    double loglk = 0.0;
    for (int64_t iPos = 0; iPos < nPos; iPos++) {
        /* Prevent underflow on large trees by comparing to maximum loglk */
        double maxloglk = -1e20;
        for (int64_t iRate = 0; iRate < (int64_t) rates.rates.size(); iRate++) {
            double site_loglk = s.site_loglk[nPos * iRate + iPos];
            if (site_loglk > maxloglk)
                maxloglk = site_loglk;
        }
        double rellk = 0; /* likelihood scaled by exp(maxloglk) */
        for (int64_t iRate = 0; iRate < (int64_t) rates.rates.size(); iRate++) {
            double lk = std::exp(s.site_loglk[nPos * iRate + iPos] - maxloglk);
            rellk += lk * dRate[iRate];
        }
        double loglk_site = maxloglk + std::log(rellk);
        loglk += loglk_site;
        if (gamma_loglk_sites != nullptr) {
            gamma_loglk_sites[iPos] = loglk_site;
        }
    }
    return loglk;
}

AbsNeighbourJoining(double)::
rescaleGammaLogLk(std::vector<numeric_t, typename op_t::Allocator> &rates, std::vector<double> &site_loglk) {
    Siteratelk s = { /*mult*/1.0, /*alpha*/1.0, rates.data(), site_loglk.data()};
    double fx, f2x;
    int i;
    fx = -gammaLogLk(s, nullptr);
    if (options.verbose > 2) {
        log << strformat("Optimizing alpha, starting at loglk %.3f", -fx) << std::endl;
    }
    for (i = 0; i < 10; i++) {
        progressReport.print("Optimizing alpha round %ld", i + 1);
        double start = fx;
        s.alpha = onedimenmin(
                0.01, s.alpha, 10.0,
                [this](double alpha, Siteratelk &s) {/* OptAlpha */
                    s.alpha = alpha;
                    return (-gammaLogLk(s, nullptr));
                },
                s, 0.001, 0.001, fx, f2x);
        if (options.verbose > 2) {
            log << strformat("Optimize alpha round %d to %.3f lk %.3f", i + 1, s.alpha, -fx) << std::endl;
        }
        s.mult = onedimenmin(
                0.01, s.mult, 10.0,
                [this](double mult, Siteratelk &s) {/* OptMult */
                    s.mult = mult;
                    return (-gammaLogLk(s, nullptr));
                },
                s, 0.001, 0.001, fx, f2x);
        if (options.verbose > 2) {
            log << strformat("Optimize mult round %d to %.3f lk %.3f", i + 1, s.mult, -fx) << std::endl;
        }
        if (fx > start - 0.001) {
            if (options.verbose > 2) {
                log << "Optimizing alpha & mult converged" << std::endl;
            }
            break;
        }
    }

    std::vector<double> gamma_loglk_sites(nPos);
    double _gammaLogLk = gammaLogLk(s, /*OUT*/gamma_loglk_sites.data());
    if (options.verbose > 0) {
        log << strformat("Gamma(%d) LogLk = %.3f alpha = %.3f rescaling lengths by %.3f",
                         options.nRateCats, _gammaLogLk, s.alpha, 1 / s.mult) << std::endl;
    }
    if (!options.logFileName.empty()) {
        log << strformat("Gamma%dLogLk\t%.3f\tApproximate\tAlpha\t%.3f\tRescale\t%.3f",
                         options.nRateCats, _gammaLogLk, s.alpha, 1 / s.mult) << std::endl;
        log << strformat("Gamma%d\tSite\tLogLk", options.nRateCats);
        for (int64_t iRate = 0; iRate < options.nRateCats; iRate++) {
            log << strformat("\tr=%.3f", rates[iRate] / s.mult);
        }
        log << std::endl;

        for (int64_t iPos = 0; iPos < nPos; iPos++) {
            log << strformat("Gamma%d\t%ld\t%.3f", options.nRateCats, iPos, gamma_loglk_sites[iPos]);
            for (int64_t iRate = 0; iRate < options.nRateCats; iRate++) {
                log << strformat("\t%.3f", site_loglk[nPos * iRate + iPos]);
            }
            log << std::endl;
        }
    }
    return (1.0 / s.mult);
}


AbsNeighbourJoining(double)::pGamma(double x, double alpha) {
    /* scale = 1/alpha */
    return incompleteGamma(x * alpha, alpha, lnGamma(alpha));
}

AbsNeighbourJoining(void)::MLSiteRates(std::vector<numeric_t, typename op_t::Allocator> &_rates) {
    _rates.resize(options.nRateCats);
    /* Even spacing from 1/nRate to nRate */
    double logNCat = std::log((double) options.nRateCats);
    double logMinRate = -logNCat;
    double logMaxRate = logNCat;
    double logd = (logMaxRate - logMinRate) / (double) (options.nRateCats - 1);

    for (int64_t i = 0; i < options.nRateCats; i++) {
        _rates[i] = std::exp(logMinRate + logd * (double) i);
    }
}


AbsNeighbourJoining(void)::MLSiteLikelihoodsByRate(std::vector<numeric_t, typename op_t::Allocator> &_rates,
                                                   std::vector<double> &site_loglk) {
    site_loglk.resize(nPos * options.nRateCats);
    /* save the original rates */
    assert(!rates.rates.empty());
    std::vector<numeric_t, typename op_t::Allocator> oldRates(rates.rates);

    /* Compute site likelihood for each rate */
    for (int64_t iRate = 0; iRate < options.nRateCats; iRate++) {
        for (int64_t i = 0; i < (int64_t) rates.rates.size(); i++) {
            rates.rates[i] = _rates[iRate];
        }
        recomputeMLProfiles();
        double loglk = treeLogLk(&site_loglk[nPos * iRate]);
        progressReport.print("Site likelihoods with rate category %ld of %d", iRate + 1, options.nRateCats);
        if (options.verbose > 2) {
            log << strformat("Rate %.3f Loglk %.3f SiteLogLk", _rates[iRate], loglk);
            for (int64_t iPos = 0; iPos < nPos; iPos++) {
                log << strformat("\t%.3f", site_loglk[nPos * iRate + iPos]);
            }
            log << std::endl;
        }
    }

    /* restore original rates and profiles */
    rates.rates = std::move(oldRates);
    recomputeMLProfiles();
}

AbsNeighbourJoining(double)::MLQuartetLogLk(Profile &pA, Profile &pB, Profile &pC, Profile &pD, double branchLengths[5],
                                            double siteLikelihoods[]) {
    Profile pAB(nPos, /*nConstraints*/0);
    Profile pCD(nPos, /*nConstraints*/0);

    posteriorProfile(pAB, pA, pB, branchLengths[0], branchLengths[1]);
    posteriorProfile(pCD, pC, pD, branchLengths[2], branchLengths[3]);
    if (siteLikelihoods != nullptr) {
        for (int64_t i = 0; i < nPos; i++) {
            siteLikelihoods[i] = 1.0;
        }
    }
    /* Roughly, P(A,B,C,D) = P(A) P(B|A) P(D|C) P(AB | CD) */
    return pairLogLk(pA, pB, branchLengths[0] + branchLengths[1], siteLikelihoods)
           + pairLogLk(pC, pD, branchLengths[2] + branchLengths[3], siteLikelihoods)
           + pairLogLk(pAB, pCD, branchLengths[4], siteLikelihoods);

}

AbsNeighbourJoining(void)::setMLRates() {
    assert(options.nRateCats > 0);
    rates.reset(1, nPos);

    if (options.nRateCats == 1) {
        recomputeMLProfiles();
        return;
    }
    std::vector<numeric_t, typename op_t::Allocator> _rates;
    std::vector<double> site_loglk;

    MLSiteRates(_rates);
    MLSiteLikelihoodsByRate(_rates, site_loglk);

    /* Select best rate for each site, correcting for the prior
       For a prior, use a gamma distribution with shape parameter 3, scale 1/3, so
       Prior(rate) ~ rate**2 * exp(-3*rate)
       log Prior(rate) = C + 2 * log(rate) - 3 * rate
    */
    double sumRates = 0;
    for (int64_t iPos = 0; iPos < nPos; iPos++) {
        int64_t iBest = -1;
        double dBest = -1e20;
        for (int64_t iRate = 0; iRate < options.nRateCats; iRate++) {
            double site_loglk_with_prior = site_loglk[nPos * iRate + iPos]
                                           + 2.0 * std::log(_rates[iRate]) - 3.0 * _rates[iRate];
            if (site_loglk_with_prior > dBest) {
                iBest = iRate;
                dBest = site_loglk_with_prior;
            }
        }
        if (options.verbose > 2) {
            log << strformat("Selected rate category %ld rate %.3f for position %ld",
                             iBest, _rates[iBest], iPos + 1) << std::endl;
        }
        rates.ratecat[iPos] = iBest;
        sumRates += _rates[iBest];
    }

    /* Force the rates to average to 1 */
    double avgRate = sumRates / nPos;
    for (int64_t iRate = 0; iRate < options.nRateCats; iRate++) {
        _rates[iRate] /= avgRate;
    }

    /* Save the rates */
    rates.rates = std::move(_rates);

    /* Update profiles based on rates */
    recomputeMLProfiles();

    if (options.verbose) {
        log << strformat("Switched to using %d rate categories (CAT approximation)", options.nRateCats) << std::endl;
        log << strformat("Rate categories were divided by %.3f so that average rate = 1.0", avgRate) << std::endl;
        log << strformat("CAT-based log-likelihoods may not be comparable across runs") << std::endl;
        if (!options.gammaLogLk) {
            log << strformat("Use -gamma for approximate but comparable Gamma(20) log-likelihoods") << std::endl;
        }
    }
}

AbsNeighbourJoining(void)::readTreeError(const std::string &err, const std::string &token) {
    throw std::invalid_argument(strformat("Tree parse error: unexpected token '%s' -- %s",
                                          token.empty() ? "(End of file)" : token.c_str(),
                                          err.c_str()
    ));
}

AbsNeighbourJoining(void)::logMLRates() {
    if (!options.logFileName.empty()) {
        log << "NCategories" << rates.rates.size() << std::endl;
        log << "Rates";

        assert(!rates.rates.empty());
        for (int64_t iRate = 0; iRate < (int64_t) rates.rates.size(); iRate++) {
            log << strformat(" %f", rates.rates[iRate]);
        }
        log << std::endl;
        log << "SiteCategories";
        for (int64_t iPos = 0; iPos < nPos; iPos++) {
            int64_t iRate = rates.ratecat[iPos];
            log << " " << iRate + 1;
        }
        log << std::endl;
    }
}

AbsNeighbourJoining(void)::logTree(const std::string &format, int64_t i, std::vector<std::string> &names,
                                   Uniquify &unique) {
    if (!options.logFileName.empty()) {
        log << strformat(format, i) << "\t";
        printNJ(log, names, unique,/*support*/false);
    }
}


AbsNeighbourJoining(void)::printVisualTree(int node, bool isFirst, const std::string &prefix) {
    log << prefix << (isFirst ? "├──" : "└──") << node << std::endl;

    if (child[node].nChild > 0) {
        printVisualTree(child[node].child[0], child[node].nChild > 1, prefix + (isFirst ? "│  " : "   "));
        if (child[node].nChild > 1) {
            printVisualTree(child[node].child[1], child[node].nChild > 2, prefix + (isFirst ? "│  " : "   "));
            if (child[node].nChild > 2) {
                printVisualTree(child[node].child[2], false, prefix + (isFirst ? "│  " : "   "));
            }
        }
    }

}

AbsNeighbourJoining(int64_t)::
treeChunks(std::vector<int64_t> &partition, std::vector<int64_t> &weights, std::vector<std::vector<int64_t>> &chunks) {
    chunks.resize(options.threads);
    std::vector<int64_t> boxWeights(options.threads, 0);
    int64_t sum = 0;
    for (int i = 0; i < (int) partition.size(); i++) {
        chunks[0].push_back(partition[i]);
        sum += weights[partition[i]];
        boxWeights[0] += weights[partition[i]];
        for (int j = 1; j < options.threads; j++) {
            if (boxWeights[j - 1] > boxWeights[j]) {
                std::swap(chunks[j - 1], chunks[j]);
                std::swap(boxWeights[j - 1], boxWeights[j]);
            } else {
                break;
            }
        }
    }
    return sum - boxWeights.back();
}

/*
 * The weight of a partition is the sum of the thread with the most nodes and the excluded nodes, that will be
 * processed sequentially.
 * */
AbsNeighbourJoining(int64_t)::treeWeight(std::vector<int64_t> &partition, std::vector<int64_t> &weights,
                                         int64_t rootWeight) {
    if (options.threadSubtrees < (int64_t) partition.size()) {
        partition.resize(options.threadSubtrees);
    }
    if (partition.empty()) {
        return rootWeight;
    } else if (partition.size() <= (size_t) options.threads) {
        int64_t weight = 0;
        for (int64_t branch: partition) {
            weight += weights[branch];
        }
        return rootWeight - weight + weights[partition[0]];
    } else {
        std::vector<std::vector<int64_t>> chunks;
        return rootWeight - treeChunks(partition, weights, chunks);
    }
}

AbsNeighbourJoining(int64_t)::treePenWeight(int64_t branch, std::vector<int64_t> &weights, int deep) {
    std::vector<int64_t> nodes, nodes2;
    nodes.push_back(branch);

    for (int i = 0; i < deep; i++) {
        for (int64_t node: nodes) {
            for (int j = 0; j < child[node].nChild; j++) {
                nodes2.push_back(child[node].child[j]);
            }
        }
        nodes.clear();
        std::swap(nodes, nodes2);
        if (nodes.empty()) {
            return 0;
        }
    }
    int64_t weight = 0;
    for (int64_t node: nodes) {
        weight += weights[node];
    }
    return weight;
}

AbsNeighbourJoining(std::vector<int64_t>)::treeChunksPaths() {
    std::vector<int64_t> path;
    for (auto &chunk: chunksPartitionCache) {
        for (int64_t branchRoot: chunk) {
            int64_t node = branchRoot;
            path.push_back(node);
            while (node != -1) {
                node = parent[node];
                path.push_back(node);
            }
        }
    }
    return path;
}

AbsNeighbourJoining(void)::treePartition(std::vector<std::vector<int64_t>> &chunks, std::vector<ibool> &traversal,
                                         int deepPen) {
    if (!chunksToRootCache.empty() && treeChunksPaths() == chunksToRootCache) {
        if (options.verbose > 0 && options.threadsVerbose) {
            log << "The tree has not changed and it was divided using the previous subtrees" << std::endl;
        }
        chunks = chunksPartitionCache;
        return;
    }
    std::vector<ibool> is_root(maxnode, false);
    std::vector<int64_t> weights(maxnode, 1);
    weights[root] = maxnode;
    std::vector<int64_t> pen_weights(maxnode, deepPen > 0 ? 0 : 1);

    std::list<int64_t> branchs;
    auto cmp = [&weights](int64_t x, int64_t y) { return weights[x] > weights[y]; };
    std::vector<int64_t> sortedPartition;
    std::vector<int64_t> partition;
    std::vector<int64_t> bestPartition;
    int64_t bestWeight = weights[root];

    for (int64_t i = 0; i < maxnode; i++) {
        if (child[i].nChild == 0) {
            branchs.push_back(i);
            sortedPartition.push_back(i);
            is_root[i] = true;
        }
    }

    while (true) {
        int64_t node = branchs.front();
        int64_t pnode = parent[node];
        branchs.pop_front();
        /* if the parent is root and the nodes are its children, we are finished */
        if (pnode == root) {
            branchs.push_back(node);
            if (branchs.size() < 4) {
                int64_t sum = 1;
                for (int64_t ni: branchs) {
                    sum += weights[ni];
                }
                if (weights[root] == sum) {
                    break;
                }
            }
            continue;
        }

        /*if it is an only child, add its father*/
        if (child[pnode].nChild == 1) {
            is_root[node] = false;
            is_root[pnode] = true;
            weights[pnode] = weights[node] + 1;
            branchs.push_back(pnode);
            continue;
        }

        int64_t sib = child[pnode].child[0] != node ? child[pnode].child[0] : child[pnode].child[1];
        /*we can only proceed if its sibling is a root*/
        if (!is_root[sib]) {
            /*if sibling has not processed the parent, it must keep its nodes.*/
            if (is_root[node]) {
                branchs.push_back(node);
            }
            continue;
        }

        is_root[sib] = false;
        is_root[node] = false;
        is_root[pnode] = true;
        weights[pnode] = weights[node] + weights[sib] + 1;
        branchs.push_back(pnode);

        auto it = std::lower_bound(sortedPartition.begin(), sortedPartition.end(), node, cmp);
        while (*it != node) {
            it++;
        }
        sortedPartition.erase(it);
        it = std::lower_bound(sortedPartition.begin(), sortedPartition.end(), sib, cmp);
        while (*it != sib) {
            it++;
        }
        sortedPartition.erase(it);
        sortedPartition.insert(std::lower_bound(sortedPartition.begin(), sortedPartition.end(), pnode, cmp), pnode);

        pen_weights[node] = pen_weights[sib] = 0; /*sibling is still in the partition, so we delete children nodes*/
        pen_weights[pnode] = treePenWeight(pnode, weights, deepPen);

        for (int64_t branch: sortedPartition) {
            if (pen_weights[branch] > 0) {
                partition.push_back(branch);
            }
        }

        int64_t pw = treeWeight(partition, pen_weights, weights[root]);
        if (pw < bestWeight && !partition.empty()) {
            bestWeight = pw;
            std::swap(bestPartition, partition);
        }
        partition.clear();
    }

    for (int64_t branch: bestPartition) {
        pen_weights[branch] = treePenWeight(branch, weights, deepPen);
    }

    treeChunks(bestPartition, pen_weights, chunks);
    chunksPartitionCache = chunks;
    chunksToRootCache = treeChunksPaths();

    if (options.verbose > 0 && options.threadsVerbose) {
        log << "The tree has " << weights[root] << " nodes and it was divided into " << bestPartition.size()
            << " subtrees:" << std::endl;
        int64_t skipped = weights[root];
        for (int i = 0; i < (int) chunks.size(); i++) {
            int64_t w = 0;
            for (int j = 0; j < (int) chunks[i].size(); j++) {
                w += pen_weights[chunks[i][j]];
            }
            skipped -= w;
            log << strformat("    thread%2d(%3.2f%%):", i, 100.0 * w / weights[root]);
            log << "branchs[";
            for (int j = 0; j < (int) chunks[i].size(); j++) {
                log << chunks[i][j];
                if (j + 1 < (int) chunks[i].size()) {
                    log << ", ";
                }
            }
            log << "], nodes " << w << std::endl;
        }
        log << strformat("    skipped (%3.2f%%): nodes %d", skipped * 100.0 / weights[root], skipped);
        log << std::endl;
        log << strformat(" total (%3.2f%%), nodes %d, theoretical speedup %.2f of %d Time %.2f",
                         100.0 * bestWeight / weights[root],
                         bestWeight, weights[root] / (float) bestWeight, options.threads, progressReport.clockDiff());
        log << std::endl;
    }

}

#ifdef __GNUC__ //disable multiple false warning in gcc 8+
    #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

AbsNeighbourJoining(inline int64_t)::traverseNNI(int64_t iRound, int64_t nRounds, int64_t &nNNIThisRound, bool useML,
                                                 std::vector<NNIStats> &stats, double &dMaxDelta, int64_t node,
                                                 std::unique_ptr<Profile> upProfiles[],
                                                 std::vector<ibool> &traversal) {
    const bool parallel = omp_in_parallel();
    double supportThreshold = useML ? Constants::treeLogLkDelta : options.MEMinDelta;
    int64_t iDone = 0;
    int64_t branchRoot = node;
    bool bUp;
    while ((node = traversePostorder(node, traversal, &bUp, branchRoot)) >= 0) {//level-2
        if (node < (int64_t) nSeqs || node == root) {
            continue; /* nothing to do for leaves or root */
        }
        if (bUp) {
            if (options.verbose > 2) {
                log << "Going up back to node " << node << std::endl;
            }
            /* No longer needed */
            for (int64_t i = 0; i < child[node].nChild; i++) {
                upProfiles[child[node].child[i]].reset();
            }
            upProfiles[node].reset();
            recomputeProfile(upProfiles, node, useML);
            continue;
        }
        if ((iDone % 100) == 0 && !parallel) {
            std::string buf;
            buf.reserve(100);
            buf += useML ? "ML" : "ME";
            buf += " NNI round %ld of %ld, %ld of %ld splits";
            if (iDone > 0) {
                buf += strformat(", %ld changes", nNNIThisRound);
            }
            if (nNNIThisRound > 0) {
                buf += strformat(" (max delta %.3f)", dMaxDelta);
            }
            progressReport.print(buf, iRound + 1, nRounds, iDone + 1, (int64_t) (maxnode - nSeqs));
        }
        iDone++;

        Profile *profiles4[4] = {nullptr, nullptr, nullptr, nullptr};
        int64_t nodeABCD[4];
        /* Note -- during the first round of ML NNIs, we use the min-evo-based branch lengths,
           which may be suboptimal */
        setupABCD(node, profiles4, upProfiles, nodeABCD, useML);

        /* Given our 4 profiles, consider doing a swap */
        int64_t nodeA = nodeABCD[0];
        int64_t nodeB = nodeABCD[1];
        int64_t nodeC = nodeABCD[2];
        int64_t nodeD = nodeABCD[3];

        auto choice = ABvsCD;

        if (options.verbose > 2) {
            log << strformat("Considering NNI around %ld: Swap A=%ld B=%ld C=%ld D=up(%ld) or parent %ld",
                             node, nodeA, nodeB, nodeC, nodeD, parent[node]) << std::endl;
        }
        if (options.verbose > 3 && useML) {
            double len[5] = {branchlength[nodeA], branchlength[nodeB], branchlength[nodeC],
                             branchlength[nodeD],
                             branchlength[node]};
            for (int i = 0; i < 5; i++) {
                if (len[i] < options.MLMinBranchLength)
                    len[i] = options.MLMinBranchLength;
            }
            log << strformat("Starting quartet likelihood %.3f len %.3f %.3f %.3f %.3f %.3f",
                             MLQuartetLogLk(*profiles4[0], *profiles4[1], *profiles4[2], *profiles4[3], len, /*site_lk*/
                                            nullptr),
                             len[0], len[1], len[2], len[3], len[4]) << std::endl;
        }

        numeric_t newlength[5];
        double criteria[3];
        if (useML) {
            for (int i = 0; i < 4; i++) {
                newlength[i] = branchlength[nodeABCD[i]];
            }
            newlength[4] = branchlength[node];
            bool bFast = options.mlAccuracy < 2 && stats[node].age > 0;
            choice = MLQuartetNNI(profiles4, criteria, newlength, bFast);
        } else {
            choice = chooseNNI(profiles4, criteria);
            /* invert criteria so that higher is better, as in ML case, to simplify code below */
            for (int i = 0; i < 3; i++)
                criteria[i] = -criteria[i];
        }

        if (choice == ACvsBD) {
            /* swap B and C */
            replaceChild(node, nodeB, nodeC);
            replaceChild(parent[node], nodeC, nodeB);
        } else if (choice == ADvsBC) {
            /* swap A and C */
            replaceChild(node, nodeA, nodeC);
            replaceChild(parent[node], nodeC, nodeA);
        }

        if (useML) {
            /* update branch length for the internal branch, and of any
           branches that lead to leaves, b/c those will not are not
           the internal branch for NNI and would not otherwise be set.
            */
            if (choice == ADvsBC) {
                /* For ADvsBC, MLQuartetNNI swaps B with D, but we swap A with C */
                double length2[5] = {newlength[LEN_C], newlength[LEN_D],
                                     newlength[LEN_A], newlength[LEN_B],
                                     newlength[LEN_I]};

                for (int i = 0; i < 5; i++) {
                    newlength[i] = length2[i];
                }
                /* and swap A and C */
                double tmp = newlength[LEN_A];
                newlength[LEN_A] = newlength[LEN_C];
                newlength[LEN_C] = tmp;
            } else if (choice == ACvsBD) {
                /* swap B and C */
                double tmp = newlength[LEN_B];
                newlength[LEN_B] = newlength[LEN_C];
                newlength[LEN_C] = tmp;
            }

            branchlength[node] = newlength[LEN_I];
            branchlength[nodeA] = newlength[LEN_A];
            branchlength[nodeB] = newlength[LEN_B];
            branchlength[nodeC] = newlength[LEN_C];
            branchlength[nodeD] = newlength[LEN_D];
        }

        if (options.verbose > 2 && (choice != ABvsCD || options.verbose > 2)) {
            log << strformat("NNI around %ld: Swap A=%ld B=%ld C=%ld D=out(C) -- choose %s %s %.4f",
                             node, nodeA, nodeB, nodeC,
                             choice == ACvsBD ? "AC|BD" : (choice == ABvsCD ? "AB|CD" : "AD|BC"),
                             useML ? "delta-loglk" : "-deltaLen",
                             criteria[choice] - criteria[ABvsCD]) << std::endl;
        }
        if (options.verbose >= 3 && options.slow && useML) {
            log << strformat("Old tree lk -- %.4f\n", treeLogLk(/*site_likelihoods*/nullptr)) << std::endl;
        }

        /* update stats, *dMaxDelta, etc. */
        if (choice == ABvsCD) {
            stats[node].age++;
        } else {
            if (useML) {
                options.debug.nML_NNI++;
            } else {
                options.debug.nNNI++;
            }
            nNNIThisRound++;
            stats[node].age = 0;
            stats[nodeA].age = 0;
            stats[nodeB].age = 0;
            stats[nodeC].age = 0;
            stats[nodeD].age = 0;
        }
        stats[node].delta = criteria[choice] - criteria[ABvsCD]; /* 0 if ABvsCD */
        if (stats[node].delta > dMaxDelta) {
            dMaxDelta = stats[node].delta;
        }

        /* support is improvement of score for self over better of alternatives */
        stats[node].support = 1e20;
        for (int i = 0; i < 3; i++) {
            if (choice != i && criteria[choice] - criteria[i] < stats[node].support) {
                stats[node].support = criteria[choice] - criteria[i];
            }
        }

        /* subtreeAge is the number of rounds since self or descendent had a significant improvement */
        if (stats[node].delta > supportThreshold) {
            stats[node].subtreeAge = 0;
        } else {
            stats[node].subtreeAge++;
            for (int i = 0; i < 2; i++) {
                int64_t ichild = child[node].child[i];
                if (stats[node].subtreeAge > stats[ichild].subtreeAge) {
                    stats[node].subtreeAge = stats[ichild].subtreeAge;
                }
            }
        }

        /* update profiles and free up unneeded up-profiles */
        if (choice == ABvsCD) {
            /* No longer needed */
            upProfiles[nodeA].reset();
            upProfiles[nodeB].reset();
            upProfiles[nodeC].reset();
            recomputeProfile(upProfiles, node, useML);
            if (options.slow && useML) {
                updateForNNI(node, upProfiles, useML);
            }
        } else {
            updateForNNI(node, upProfiles, useML);
        }
        if (options.verbose > 2 && options.slow && useML) {
            /* Note we recomputed profiles back up to root already if slow */
            printNJInternal(log, /*useLen*/true);
            log << strformat("New tree lk -- %.4f\n", treeLogLk(/*site_likelihoods*/nullptr));
        }
    }
    return iDone;
}

AbsNeighbourJoining(int64_t)::DoNNI(int64_t iRound, int64_t nRounds, bool useML, std::vector<NNIStats> &stats,
                                    double &dMaxDelta) {
/* For each non-root node N, with children A,B, sibling C, and uncle D,
     we compare the current topology AB|CD to the alternate topologies
     AC|BD and AD|BC, by using the 4 relevant profiles.

     If useML is true, it uses quartet maximum likelihood, and it
     updates branch lengths as it goes.

     If useML is false, it uses the minimum-evolution criterion with
     log-corrected distances on profiles.  (If logdist is false, then
     the log correction is not done.) If useML is false, then NNI()
     does NOT modify the branch lengths.

     Regardless of whether it changes the topology, it recomputes the
     profile for the node, using the pairwise distances and BIONJ-like
     weightings (if bionj is set). The parent's profile has changed,
     but recomputing it is not necessary because we will visit it
     before we need it (we use postorder, so we may visit the sibling
     and its children before we visit the parent, but we never
     consider an ancestor's profile, so that is OK). When we change
     the parent's profile, this alters the uncle's up-profile, so we
     remove that.  Finally, if the topology has changed, we remove the
     up-profiles of the nodes.

     If we do an NNI during post-order traversal, the result is a bit
     tricky. E.g. if we are at node N, and have visited its children A
     and B but not its uncle C, and we do an NNI that swaps B & C,
     then the post-order traversal will visit C, and its children, but
     then on the way back up, it will skip N, as it has already
     visited it.  So, the profile of N will not be recomputed: any
     changes beneath C will not be reflected in the profile of N, and
     the profile of N will be slightly stale. This will be corrected
     on the next round of NNIs.
  */
    double supportThreshold = useML ? Constants::treeLogLkDelta : options.MEMinDelta;
    int64_t nNNIThisRound = 0;
    dMaxDelta = 0.0;

    if (nSeqs <= 3) {
        return 0;            /* nothing to do */
    }
    if (options.verbose > 2) {
        log << strformat("Beginning round %ld of NNIs with ml? %d", iRound, useML ? 1 : 0) << std::endl;
        printNJInternal(log, /*useLen*/useML && iRound > 0 ? true : false);
    }

    /* For each node the upProfile or NULL */
    std::vector<std::unique_ptr<Profile>> upProfiles(maxnodes);
    std::vector<ibool> traversal(maxnodes, false);


    /* Identify nodes we can skip traversing into */
    if (options.fastNNI) {
        for (int64_t node = 0; node < maxnode; node++) {
            if (node != root
                && node >= (int64_t) nSeqs
                && stats[node].age >= 2
                && stats[node].subtreeAge >= 2
                && stats[node].support > supportThreshold) {
                int64_t nodeABCD[4];
                setupABCD(node, nullptr, nullptr, nodeABCD, useML);
                int i;
                for (i = 0; i < 4; i++) {
                    if (stats[nodeABCD[i]].age == 0 && stats[nodeABCD[i]].support > supportThreshold) {
                        break;
                    }
                }
                if (i == 4) {
                    traversal[node] = true;
                    if (options.verbose > 2) {
                        log << strformat(
                                "Skipping subtree at %ld: child %ld %ld parent %ld age %ld subtreeAge %ld support %.3f",
                                node, nodeABCD[0], nodeABCD[1], parent[node],
                                stats[node].age, stats[node].subtreeAge, stats[node].support) << std::endl;
                    }
                }
            }
        }
    }

    std::string buf = useML ? "ML" : "ME";
    buf += " NNI round %ld of %ld, %ld splits";
    progressReport.print(buf, iRound + 1, nRounds, (int64_t) (maxnode - nSeqs));

    if (options.threads > 1 && options.threadsLevel > 1) {
        std::vector<std::vector<int64_t>> chunks;
        treePartition(chunks, traversal, 2);
        int64_t done = 0;

        #pragma omp parallel
        {
            std::vector<std::unique_ptr<Profile>> upProfiles2(maxnodes);
            double dMaxDelta2 = 0.0;
            int64_t nNNIThisRound2 = 0;
            int64_t done2;

            #pragma omp for schedule(static, 1)
            for (int i = 0; i < (int) chunks.size(); i++) {
                for (int64_t subtreeRoot: chunks[i]) {
                    /*
                     * the neighborhood of a node includes the parent, so in order not to modify the subtreeRoot node,
                     * we do not process their direct children.
                     */
                    for (int j = 0; j < child[subtreeRoot].nChild; j++) {
                        int64_t directChildren = child[subtreeRoot].child[j];
                        for (int k = 0; k < child[directChildren].nChild; k++) {
                            int64_t children = child[directChildren].child[k];
                            done2 = traverseNNI(iRound, nRounds, nNNIThisRound2, useML, stats, dMaxDelta2, children,
                                                upProfiles2.data(), traversal);
                            if (options.verbose > 0) {
                                #pragma omp critical
                                {
                                    done += done2;
                                    std::string buf;
                                    buf.reserve(100);
                                    buf += useML ? "ML" : "ME";
                                    buf += " NNI round %ld of %ld, %ld of %ld splits";
                                    if (done > 0) {
                                        buf += strformat(", %ld changes", nNNIThisRound);
                                    }
                                    if (nNNIThisRound > 0) {
                                        buf += strformat(" (max delta %.3f)", dMaxDelta);
                                    }
                                    progressReport.print(buf, iRound + 1, nRounds, done + 1,
                                                         (int64_t) (maxnode - nSeqs));
                                }
                            }
                        }
                    }
                }
            }

            #pragma omp critical
            {
                if (dMaxDelta2 > dMaxDelta) {
                    dMaxDelta = dMaxDelta2;
                }
                nNNIThisRound += nNNIThisRound2;
            }
        }
    }

    traverseNNI(iRound, nRounds, nNNIThisRound, useML, stats, dMaxDelta, root, upProfiles.data(), traversal);

    if (options.verbose >= 2) {
        int nUp = 0;
        for (int64_t i = 0; i < maxnodes; i++) {
            if (upProfiles[i]) {
                nUp++;
            }
        }
        log << "N up profiles at end of NNI:  " << nUp << std::endl;
    }
    return nNNIThisRound;
}

AbsNeighbourJoining(inline int64_t)::traverseSPR(int64_t iRound, int64_t nRounds, std::unique_ptr<Profile> upProfiles[],
                                                 std::vector<ibool> &traversal, int64_t branchRoot,
                                                 double last_tot_len) {
    const bool parallel = omp_in_parallel();
    std::vector<int64_t> nodeList(maxnodes);
    int64_t nodeListLen = 0;
    int64_t node = branchRoot;
    while ((node = traversePostorder(node, traversal, /*pUp*/nullptr, branchRoot)) >= 0) {//level-4
        nodeList[nodeListLen++] = node;
    }

    assert(nodeListLen == maxnode || options.threadsLevel > 3);
    std::vector<SprStep> steps(options.maxSPRLength); /* current chain of SPRs */

    for (int64_t i = 0; i < nodeListLen; i++) {
        int64_t node = nodeList[i];
        if ((i % 100) == 0 && !parallel) {
            progressReport.print("SPR round %3ld of %3ld, %ld of %ld nodes",
                                 iRound + 1, nRounds, i + 1, nodeListLen);
        }
        if (node == root) {
            continue; /* nothing to do for root */
        }
        /* The nodes to NNI around */
        int64_t nodeAround[2] = {parent[node], sibling(node)};
        if (parent[node] == root) {
            /* NNI around both siblings instead */
            rootSiblings(node, /*OUT*/nodeAround);
        }
        bool bChanged = false;
        for (int64_t iAround = 0; iAround < 2 && bChanged == false; iAround++) {
            for (int ACFirst = 0; ACFirst < 2 && bChanged == false; ACFirst++) {
                if (options.verbose > 3) {
                    printNJInternal(log, /*useLen*/false);
                }
                int64_t chainLength = findSPRSteps(node, nodeAround[iAround], upProfiles, steps.data(), (bool) ACFirst);
                double dMinDelta = 0.0;
                int64_t iCBest = -1;
                double dTotDelta = 0.0;
                for (int64_t iC = 0; iC < chainLength; iC++) {
                    dTotDelta += steps[iC].deltaLength;
                    if (dTotDelta < dMinDelta) {
                        dMinDelta = dTotDelta;
                        iCBest = iC;
                    }
                }

                if (options.verbose > 3) {
                    log << strformat("SPR %s %ld around %ld chainLength %ld of %ld deltaLength %.5f swaps:",
                                     iCBest >= 0 ? "move" : "abandoned",
                                     node, nodeAround[iAround], iCBest + 1, chainLength, dMinDelta);
                    for (int64_t iC = 0; iC < chainLength; iC++) {
                        log << strformat(" (%ld,%ld)%.4f", steps[iC].nodes[0], steps[iC].nodes[1],
                                         steps[iC].deltaLength);
                    }
                    log << std::endl;
                }
                for (int64_t iC = chainLength - 1; iC > iCBest; iC--) {
                    unwindSPRStep(steps[iC], upProfiles);
                }
                if (options.verbose > 3) {
                    printNJInternal(log, /*useLen*/false);
                }
                while (options.slow && iCBest >= 0) {
                    double expected_tot_len = last_tot_len + dMinDelta;
                    double new_tot_len = treeLength(/*recompute*/true);
                    if (options.verbose > 2) {
                        log << strformat("Total branch-length is now %.4f was %.4f expected %.4f",
                                         new_tot_len, last_tot_len, expected_tot_len) << std::endl;
                    }
                    if (new_tot_len < last_tot_len) {
                        last_tot_len = new_tot_len;
                        break;        /* no rewinding necessary */
                    }
                    if (options.verbose > 2) {
                        log << strformat("Rewinding SPR to %ld", iCBest) << std::endl;
                    }
                    unwindSPRStep(steps[iCBest], upProfiles);
                    dMinDelta -= steps[iCBest].deltaLength;
                    iCBest--;
                }
                if (iCBest >= 0) {
                    bChanged = true;
                }
            }    /* loop over which step to take at 1st NNI */
        } /* loop over which node to pivot around */

        if (bChanged) {
            options.debug.nSPR++;        /* the SPR move is OK */
            /* make sure all the profiles are OK */
            for (int64_t j = 0; j < maxnodes; j++) {
                upProfiles[j].reset();
            }
            for (int64_t ancestor = parent[node]; ancestor >= 0; ancestor = parent[ancestor]) {
                recomputeProfile(upProfiles, ancestor, /*useML*/false);
                if (ancestor == branchRoot) {
                    break;
                }
            }
        }
    } /* end loop over subtrees to prune & regraft */
    return nodeListLen;
}

AbsNeighbourJoining(void)::SPR(int64_t iRound, int64_t nRounds) {
    /* Given a non-root node N with children A,B, sibling C, and uncle D,
       we can try to move A by doing three types of moves (4 choices):
       "down" -- swap A with a child of B (if B is not a leaf) [2 choices]
       "over" -- swap B with C
       "up" -- swap A with D
       We follow down moves with down moves, over moves with down moves, and
       up moves with either up or over moves. (Other choices are just backing
       up and hence useless.)

       As with NNIs, we keep track of up-profiles as we go. However, some of the regular
       profiles may also become "stale" so it is a bit trickier.

       We store the traversal before we do SPRs to avoid any possible infinite loop
    */
    double last_tot_len = 0.0;
    if (nSeqs <= 3 || options.maxSPRLength < 1) {
        return;
    }
    if (options.slow) {
        last_tot_len = treeLength(/*recomputeLengths*/true);
    }

    std::vector<ibool> traversal(maxnodes, false);
    std::vector<std::unique_ptr<Profile>> upProfiles(maxnodes);

    if (options.threads > 1 && options.threadsLevel > 3) {

        std::vector<std::vector<int64_t>> chunks;
        std::vector<int64_t> branchRoots;
        treePartition(chunks, traversal, options.maxSPRLength + 1);
        int64_t done = 0;

        if (options.slow) {
            lockedNodes.resize(maxnodes, false);
            for (int i = 0; i < (int) chunks.size(); i++) {
                for (int64_t subtreeRoot: chunks[i]) {
                    lockedNodes[subtreeRoot] = true;
                }
            }
        }

        #pragma omp parallel
        {
            std::vector<std::unique_ptr<Profile>> upProfiles2(maxnodes);
            int64_t done2;

            #pragma omp for schedule(static, 1)
            for (int i = 0; i < (int) chunks.size(); i++) {
                std::vector<int64_t> stackA = chunks[i];
                std::vector<int64_t> stackB;

                for (int i = 0; i < options.maxSPRLength + 1; i++) {
                    for (int64_t subtreeRoot: stackA) {
                        for (int j = 0; j < child[subtreeRoot].nChild; j++) {
                            stackB.push_back(child[subtreeRoot].child[j]);
                        }
                    }
                    stackA = std::move(stackB);
                    stackB.clear();
                    if (stackA.empty()) {
                        break;
                    }
                }

                for (int64_t subtreeRoot: stackA) {
                    done2 = traverseSPR(iRound, nRounds, upProfiles2.data(), traversal, subtreeRoot, last_tot_len);
                    if (options.verbose > 0) {
                        #pragma omp critical
                        {
                            done += done2;
                            progressReport.print("SPR round %3ld of %3ld, %ld of %ld nodes",
                                                 iRound + 1, nRounds, done + 1, maxnode);
                        }
                    }
                }

                #pragma omp critical
                {
                    branchRoots.insert(branchRoots.end(), stackA.begin(), stackA.end());
                }
            }
        }

        for (int64_t subtreeRoot: branchRoots) {
            for (int64_t j = 0; j < maxnodes; j++) {
                upProfiles[j].reset();
            }
            for (int64_t ancestor = parent[subtreeRoot]; ancestor >= 0; ancestor = parent[ancestor]) {
                recomputeProfile(upProfiles.data(), ancestor, /*useML*/false);
            }
        }
        lockedNodes.resize(0);
    }

    traverseSPR(iRound, nRounds, upProfiles.data(), traversal, root, last_tot_len);
}


AbsNeighbourJoining(double)::GTRNegLogLk(double x, GtrOpt &gtr) {
    assert(options.nCodes == 4);
    assert(gtr.iRate >= 0 && gtr.iRate < 6);
    assert(x > 0);

    double _rates[6];
    for (int i = 0; i < 6; i++) {
        _rates[i] = gtr.rates[i];
    }
    _rates[gtr.iRate] = x;

    if (options.verbose > 2) {
        log << strformat("GTR_Opt\tfreq %.5f %.5f %.5f %.5f rates %.5f %.5f %.5f %.5f %.5f %.5f",
                         gtr.freq[0], gtr.freq[1], gtr.freq[2], gtr.freq[3],
                         _rates[0], _rates[1], _rates[2], _rates[3], _rates[4], _rates[5]) << std::endl;
    }

    transmat.createGTR(options, _rates, gtr.freq); /*transmat is restored by GTRNegLogLk caller */
    recomputeMLProfiles();
    double loglk = treeLogLk(/*site_loglk*/nullptr);

    /* Do not recompute profiles -- assume the caller will do that */
    if (options.verbose > 2) {
        log << strformat("GTR LogLk(%.5f %.5f %.5f %.5f %.5f %.5f) = %f",
                         _rates[0], _rates[1], _rates[2], _rates[3], _rates[4], _rates[5], loglk) << std::endl;
    }
    return (-loglk);
}

AbsNeighbourJoining(void)::setMLGtr(double freq_in[]) {
    assert(options.nCodes == 4);
    GtrOpt gtr;
    if (freq_in != NULL) {
        for (int i = 0; i < 4; i++)
            gtr.freq[i] = freq_in[i];
    } else {
        /* n[] and sum were int in FastTree 2.1.9 and earlier -- this
           caused gtr analyses to fail on analyses with >2e9 positions */
        int64_t n[4] = {1, 1, 1, 1};    /* pseudocounts */
        for (int64_t i = 0; i < (int64_t) nSeqs; i++) {
            auto &codes = profiles[i].codes;
            for (int64_t iPos = 0; iPos < nPos; iPos++)
                if (codes[iPos] < 4) {
                    n[(int) codes[iPos]]++;
                }
        }
        int64_t sum = n[0] + n[1] + n[2] + n[3];
        for (int i = 0; i < 4; i++) {
            gtr.freq[i] = n[i] / (double) sum;
        }
    }
    for (int i = 0; i < 6; i++) {
        gtr.rates[i] = 1.0;
    }
    int64_t nRounds = options.mlAccuracy < 2 ? 2 : options.mlAccuracy;

    /*GTRNegLogLk GTRNegLogLk change the transition matrix, to avoid allocate new memory and swap them multiple times,
     * the transition matrix is copied and restore at the end */
    auto transmatCpy = transmat;

    for (int64_t i = 0; i < nRounds; i++) {
        for (gtr.iRate = 0; gtr.iRate < 6; gtr.iRate++) {
            progressReport.print("Optimizing GTR model, step %ld of %d", i * 6 + gtr.iRate + 1, 12);
            double negloglk, f2x;
            gtr.rates[gtr.iRate] = onedimenmin(/*xmin*/0.05,
                    /*xguess*/gtr.rates[gtr.iRate],
                    /*xmax*/20.0,
                    /*func*/[this](double x, GtrOpt &qo) { return GTRNegLogLk(x, qo); },
                    /*data*/gtr,
                    /*ftol*/0.001,
                    /*atol*/0.0001,
                    /*OUT*/negloglk,
                    /*OUT*/f2x);
        }
    }
    transmat = transmatCpy;

    /* normalize gtr so last rate is 1 -- specifying that rate separately is useful for optimization only */
    for (int i = 0; i < 5; i++) {
        gtr.rates[i] /= gtr.rates[5];
    }
    gtr.rates[5] = 1.0;
    if (options.verbose) {
        log << strformat("GTR Frequencies: %.4f %.4f %.4f %.4f", gtr.freq[0], gtr.freq[1], gtr.freq[2], gtr.freq[3])
            << std::endl;
        log << strformat("GTR rates(ac ag at cg ct gt) %.4f %.4f %.4f %.4f %.4f %.4f",
                         gtr.rates[0], gtr.rates[1], gtr.rates[2], gtr.rates[3], gtr.rates[4], gtr.rates[5])
            << std::endl;
    }
    transmat.createGTR(options, gtr.rates, gtr.freq);
    recomputeMLProfiles();
    optimizeAllBranchLengths();
}

/* Recomputes all branch lengths

   For internal branches such as (A,B) vs. (C,D), uses the formula

   length(AB|CD) = (d(A,C)+d(A,D)+d(B,C)+d(B,D))/4 - d(A,B)/2 - d(C,D)/2

   (where all distances are profile distances - diameters).

   For external branches (e.g. to leaves) A vs. (B,C), use the formula

   length(A|BC) = (d(A,B)+d(A,C)-d(B,C))/2
*/
AbsNeighbourJoining(inline void)::traverseUpdateBranchLengths(int64_t node, std::unique_ptr<Profile> upProfiles[],
                                                              std::vector<ibool> &traversal) {
    int64_t branchRoot = node;
    while ((node = traversePostorder(node, traversal, /*pUp*/nullptr, branchRoot)) >= 0) {//level-1
        /* reset branch length of node (distance to its parent) */
        if (node == root) {
            continue; /* no branch length to set */
        }
        if (node < (int64_t) nSeqs) { /* a leaf */
            Profile *profileA = &profiles[node];
            Profile *profileB = nullptr;
            Profile *profileC = nullptr;

            int64_t sib = sibling(node);
            if (sib == -1) { /* at root, have 2 siblings */
                int64_t sibs[2];
                rootSiblings(node, /*OUT*/sibs);
                profileB = &profiles[sibs[0]];
                profileC = &profiles[sibs[1]];
            } else {
                profileB = &profiles[sib];
                profileC = getUpProfile(upProfiles, parent[node], /*useML*/false);
            }
            Profile *profiles3[3] = {profileA, profileB, profileC};
            double d[3]; /*AB,AC,BC*/
            correctedPairDistances(profiles3, 3, /*OUT*/d);
            /* d(A,BC) = (dAB+dAC-dBC)/2 */
            branchlength[node] = (d[0] + d[1] - d[2]) / 2.0;
        } else {
            Profile *profiles4[4];
            int64_t nodeABCD[4];
            setupABCD(node, profiles4, upProfiles, nodeABCD, false);
            double d[6];
            correctedPairDistances(profiles4, 4, d);
            branchlength[node] = (d[qAC] + d[qAD] + d[qBC] + d[qBD]) / 4.0 - (d[qAB] + d[qCD]) / 2.0;

            /* no longer needed */
            upProfiles[nodeABCD[0]].reset();
            upProfiles[nodeABCD[1]].reset();
        }
    }
}


AbsNeighbourJoining(void)::updateBranchLengths() {
    if (nSeqs < 2) {
        return;
    } else if (nSeqs == 2) {
        int64_t nodeA = child[root].child[0];
        int64_t nodeB = child[root].child[1];
        Besthit h;
        profileDist(profiles[nodeA], profiles[nodeB], h);
        if (options.logdist) {
            h.dist = logCorrect(h.dist);
        }
        branchlength[nodeA] = h.dist / 2.0;
        branchlength[nodeB] = h.dist / 2.0;
        return;
    }

    std::vector<ibool> traversal(maxnodes, false);

    if (options.threads > 1 && options.threadsLevel > 0) {
        std::vector<std::vector<int64_t>> chunks;
        treePartition(chunks, traversal, 0);

        #pragma omp parallel
        {
            std::vector<std::unique_ptr<Profile>> upProfiles(maxnodes);
            #pragma omp for schedule(static, 1)
            for (int i = 0; i < (int) chunks.size(); i++) {
                for (int64_t subtreeRoot: chunks[i]) {
                    traverseUpdateBranchLengths(subtreeRoot, upProfiles.data(), traversal);
                }
            }
        }
    }

    std::vector<std::unique_ptr<Profile>> upProfiles(maxnodes);
    traverseUpdateBranchLengths(root, upProfiles.data(), traversal);
}


AbsNeighbourJoining(double)::treeLength(bool recomputeProfiles) {
    if (recomputeProfiles) {
        std::vector<ibool> traversal2(maxnodes, false);
        int64_t j = root;
        while ((j = traversePostorder(j, traversal2, /*pUp*/nullptr, root)) >= 0) {
            /* nothing to do for leaves or root */
            if (j >= (int64_t) nSeqs && j != root)
                setProfile(j, /*noweight*/-1.0);
        }
    }
    updateBranchLengths();
    double total_len = 0;
    for (int64_t iNode = 0; iNode < maxnode; iNode++)
        total_len += branchlength[iNode];
    return total_len;
}

AbsNeighbourJoining(inline void)::traverseTestSplitsMinEvo(int64_t node, SplitCount &splitcount,
                                                           std::unique_ptr<Profile> upProfiles[],
                                                           std::vector<ibool> &traversal) {
    const double tolerance = 1e-6;
    int64_t branchRoot = node;

    while ((node = traversePostorder(node, /*IN/OUT*/traversal, /*pUp*/nullptr, branchRoot)) >= 0) {//level-1
        if (node < (int64_t) nSeqs || node == root) {
            continue; /* nothing to do for leaves or root */
        }

        Profile *profiles4[4];
        int64_t nodeABCD[4];
        setupABCD(node, /*OUT*/profiles4, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, /*useML*/false);

        if (options.verbose > 2) {
            log << strformat("Testing Split around %ld: A=%ld B=%ld C=%ld D=up(%ld) or node parent %ld",
                             node, nodeABCD[0], nodeABCD[1], nodeABCD[2], nodeABCD[3], parent[node]) << std::endl;
        }

        double d[6];        /* distances, perhaps log-corrected distances, no constraint penalties */
        correctedPairDistances(profiles4, 4, /*OUT*/d);

        /* alignment-based scores for each split (lower is better) */
        double sABvsCD = d[qAB] + d[qCD];
        double sACvsBD = d[qAC] + d[qBD];
        double sADvsBC = d[qAD] + d[qBC];

        /* constraint penalties, indexed by nni_t (lower is better) */
        double p[3];
        quartetConstraintPenalties(profiles4, /*OUT*/p);

        int64_t nConstraintsViolated = 0;
        for (int64_t iC = 0; iC < (int64_t) constraintSeqs.size(); iC++) {
            if (splitViolatesConstraint(profiles4, iC)) {
                nConstraintsViolated++;
                if (options.verbose > 2) {
                    double penalty[3] = {0.0, 0.0, 0.0};
                    quartetConstraintPenaltiesPiece(profiles4, iC, /*OUT*/penalty);
                    log << strformat(
                            "Violate constraint %ld at %ld (children %ld %ld) penalties"
                            " %.3f %.3f %.3f %ld/%ld %ld/%ld %ld/%ld %ld/%ld",
                            iC, node, child[node].child[0], child[node].child[1],
                            penalty[ABvsCD], penalty[ACvsBD], penalty[ADvsBC],
                            profiles4[0]->nOn[iC], profiles4[0]->nOff[iC],
                            profiles4[1]->nOn[iC], profiles4[1]->nOff[iC],
                            profiles4[2]->nOn[iC], profiles4[2]->nOff[iC],
                            profiles4[3]->nOn[iC], profiles4[3]->nOff[iC]) << std::endl;
                }
            }
        }

        double delta = sABvsCD - std::min(sACvsBD, sADvsBC);
        bool bBadDist = delta > tolerance;
        bool bBadConstr = p[ABvsCD] > p[ACvsBD] + tolerance || p[ABvsCD] > p[ADvsBC] + tolerance;

        splitcount.nSplits++;
        if (bBadDist) {
            auto choice = sACvsBD < sADvsBC ? ACvsBD : ADvsBC;
            /* If ABvsCD is favored over the shorter NNI by constraints,
           then this is probably a bad split because of the constraint */
            if (p[choice] > p[ABvsCD] + tolerance)
                splitcount.dWorstDeltaConstrained = std::max(delta, splitcount.dWorstDeltaConstrained);
            else
                splitcount.dWorstDeltaUnconstrained = std::max(delta, splitcount.dWorstDeltaUnconstrained);
        }

        if (nConstraintsViolated > 0)
            splitcount.nConstraintViolations++; /* count splits with any violations, not #constraints in a splits */
        if (bBadDist)
            splitcount.nBadSplits++;
        if (bBadDist && bBadConstr)
            splitcount.nBadBoth++;
        if (bBadConstr && options.verbose > 2) {
            /* Which NNI would be better */
            double dist_advantage = 0;
            double constraint_penalty = 0;
            if (p[ACvsBD] < p[ADvsBC]) {
                dist_advantage = sACvsBD - sABvsCD;
                constraint_penalty = p[ABvsCD] - p[ACvsBD];
            } else {
                dist_advantage = sADvsBC - sABvsCD;
                constraint_penalty = p[ABvsCD] - p[ADvsBC];
            }
            log << strformat(
                    "Violate constraints %ld distance_advantage %.3f constraint_penalty %.3f (children %ld %ld):",
                    node, dist_advantage, constraint_penalty,
                    child[node].child[0], child[node].child[1]);
            /* list the constraints with a penalty, meaning that ABCD all have non-zero
               values and that AB|CD worse than others */
            for (int64_t iC = 0; iC < (int64_t) constraintSeqs.size(); iC++) {
                double ppart[6];
                if (quartetConstraintPenaltiesPiece(profiles4, iC, /*OUT*/ppart)) {
                    if (ppart[qAB] + ppart[qCD] > ppart[qAD] + ppart[qBC] + tolerance
                        || ppart[qAB] + ppart[qCD] > ppart[qAC] + ppart[qBD] + tolerance) {
                        log << strformat(" %ld (%ld/%ld %ld/%ld %ld/%ld %ld/%ld)", iC,
                                         profiles4[0]->nOn[iC], profiles4[0]->nOff[iC],
                                         profiles4[1]->nOn[iC], profiles4[1]->nOff[iC],
                                         profiles4[2]->nOn[iC], profiles4[2]->nOff[iC],
                                         profiles4[3]->nOn[iC], profiles4[3]->nOff[iC]);
                    }
                }
            }
            log << std::endl;
        }
        /* no longer needed */
        upProfiles[nodeABCD[0]].reset();
        upProfiles[nodeABCD[1]].reset();
    }
}

AbsNeighbourJoining(void)::testSplitsMinEvo(SplitCount &splitcount) {
    splitcount.nBadSplits = 0;
    splitcount.nConstraintViolations = 0;
    splitcount.nBadBoth = 0;
    splitcount.nSplits = 0;
    splitcount.dWorstDeltaUnconstrained = 0.0;
    splitcount.dWorstDeltaConstrained = 0.0;

    std::vector<ibool> traversal(maxnodes, false);

    if (options.threads > 1 && options.threadsLevel > 0) {
        std::vector<std::vector<int64_t>> chunks;
        treePartition(chunks, traversal, 0);

        #pragma omp parallel
        {
            std::vector<std::unique_ptr<Profile>> upProfiles2(maxnodes);
            SplitCount splitcount2 = splitcount;

            #pragma omp for schedule(static, 1)
            for (int i = 0; i < (int) chunks.size(); i++) {
                for (int64_t subtreeRoot: chunks[i]) {
                    traverseTestSplitsMinEvo(subtreeRoot, splitcount2, upProfiles2.data(), traversal);
                }
            }
            #pragma omp critical
            {
                splitcount.nSplits += splitcount2.nSplits;
                splitcount.nConstraintViolations += splitcount2.nSplits;
                splitcount.nBadSplits += splitcount2.nSplits;
                splitcount.nBadBoth += splitcount2.nSplits;
                splitcount.dWorstDeltaConstrained = std::max(splitcount2.dWorstDeltaConstrained,
                                                             splitcount.dWorstDeltaConstrained);
                splitcount.dWorstDeltaUnconstrained = std::max(splitcount2.dWorstDeltaUnconstrained,
                                                               splitcount.dWorstDeltaUnconstrained);
            }
        }
    }

    std::vector<std::unique_ptr<Profile>> upProfiles(maxnodes);
    traverseTestSplitsMinEvo(root, splitcount, upProfiles.data(), traversal);

}

AbsNeighbourJoining(void)::testSplitsML(SplitCount &splitcount) {
    splitcount.nBadSplits = 0;
    splitcount.nConstraintViolations = 0;
    splitcount.nBadBoth = 0;
    splitcount.nSplits = 0;
    splitcount.dWorstDeltaUnconstrained = 0;
    splitcount.dWorstDeltaConstrained = 0;


    std::vector<ibool> traversal(maxnodes, false);

    std::vector<int64_t> col;
    if (options.nBootstrap > 0) {
        resampleColumns(col);
    }

    if (options.threads > 1 && options.threadsLevel > 0) {
        std::vector<std::vector<int64_t>> chunks;
        treePartition(chunks, traversal, 0);
        int64_t done = 0;

        #pragma omp parallel
        {
            std::vector<std::unique_ptr<Profile>> upProfiles2(maxnodes);
            SplitCount splitcount2 = splitcount;
            int64_t done2;

            #pragma omp for schedule(static, 1)
            for (int i = 0; i < (int) chunks.size(); i++) {
                for (int64_t subtreeRoot: chunks[i]) {
                    done2 = traverseTestSplitsML(subtreeRoot, splitcount2, col, upProfiles2.data(), traversal);
                    if (options.verbose > 0) {
                        #pragma omp critical
                        {
                            done += done2;
                            progressReport.print("ML split tests for %6ld of %6ld internal splits", done,
                                                 (int64_t) (nSeqs - 3));
                        }
                    }
                }
            }
            #pragma omp critical
            {
                splitcount.nSplits += splitcount2.nSplits;
                splitcount.nConstraintViolations += splitcount2.nSplits;
                splitcount.nBadSplits += splitcount2.nSplits;
                splitcount.nBadBoth += splitcount2.nSplits;
                splitcount.dWorstDeltaConstrained = std::max(splitcount2.dWorstDeltaConstrained,
                                                             splitcount.dWorstDeltaConstrained);
                splitcount.dWorstDeltaUnconstrained = std::max(splitcount2.dWorstDeltaUnconstrained,
                                                               splitcount.dWorstDeltaUnconstrained);
            }
        }
    }

    std::vector<std::unique_ptr<Profile>> upProfiles(maxnodes);
    traverseTestSplitsML(root, splitcount, col, upProfiles.data(), traversal);
}

AbsNeighbourJoining(inline int64_t)::traverseTestSplitsML(int64_t node, SplitCount &splitcount,
                                                          const std::vector<int64_t> &col,
                                                          std::unique_ptr<Profile> upProfiles[],
                                                          std::vector<ibool> &traversal) {
    const bool parallel = omp_in_parallel();
    const double tolerance = 1e-6;
    std::vector<double> site_likelihoods(3 * nPos);

    int64_t iNodesDone = 0;
    while ((node = traversePostorder(node, /*IN/OUT*/traversal, /*pUp*/nullptr, root)) >= 0) {//level-1
        if (node < (int64_t) nSeqs || node == root) {
            continue; /* nothing to do for leaves or root */
        }

        if (iNodesDone > 0 && (iNodesDone % 100) == 0 && !parallel) {
            progressReport.print("ML split tests for %6ld of %6ld internal splits", iNodesDone, (int64_t) (nSeqs - 3));
        }
        iNodesDone++;

        Profile *profiles4[4];
        int64_t nodeABCD[4];
        setupABCD(node, /*OUT*/profiles4, /*IN/OUT*/upProfiles, /*OUT*/nodeABCD, /*useML*/true);
        double loglk[3];
        double len[5];
        for (int i = 0; i < 4; i++) {
            len[i] = branchlength[nodeABCD[i]];
        }
        len[4] = branchlength[node];
        double lenABvsCD[5] = {len[LEN_A], len[LEN_B], len[LEN_C], len[LEN_D], len[LEN_I]};
        double lenACvsBD[5] = {len[LEN_A], len[LEN_C], len[LEN_B], len[LEN_D], len[LEN_I]};   /* Swap B & C */
        double lenADvsBC[5] = {len[LEN_A], len[LEN_D], len[LEN_C], len[LEN_B], len[LEN_I]};   /* Swap B & D */

        if (!parallel) {
            #pragma omp parallel
            #pragma omp sections
            {
                #pragma omp section
                {
                    /* Lengths are already optimized for ABvsCD */
                    loglk[ABvsCD] = MLQuartetLogLk(*profiles4[0], *profiles4[1], *profiles4[2], *profiles4[3],
                                                   lenABvsCD, &site_likelihoods[ABvsCD * nPos]);
                }
                #pragma omp section
                {
                    loglk[ACvsBD] = MLQuartetOptimize(*profiles4[0], *profiles4[2], *profiles4[1], *profiles4[3],
                                                      lenACvsBD, /*pStarTest*/nullptr,
                                                      &site_likelihoods[ACvsBD * nPos]);
                }
                #pragma omp section
                {
                    loglk[ADvsBC] = MLQuartetOptimize(*profiles4[0], *profiles4[3], *profiles4[2], *profiles4[1],
                                                      lenADvsBC, /*pStarTest*/nullptr,
                                                      &site_likelihoods[ADvsBC * nPos]);
                }
            }
        } else {
            /* Lengths are already optimized for ABvsCD */
            loglk[ABvsCD] = MLQuartetLogLk(*profiles4[0], *profiles4[1], *profiles4[2], *profiles4[3],
                                           lenABvsCD, &site_likelihoods[ABvsCD * nPos]);

            loglk[ACvsBD] = MLQuartetOptimize(*profiles4[0], *profiles4[2], *profiles4[1], *profiles4[3],
                                              lenACvsBD, /*pStarTest*/nullptr,
                                              &site_likelihoods[ACvsBD * nPos]);

            loglk[ADvsBC] = MLQuartetOptimize(*profiles4[0], *profiles4[3], *profiles4[2], *profiles4[1],
                                              lenADvsBC, /*pStarTest*/nullptr,
                                              &site_likelihoods[ADvsBC * nPos]);
        }

        /* do a second pass on the better alternative if it is close */
        if (loglk[ACvsBD] > loglk[ADvsBC]) {
            if (options.mlAccuracy > 1 || loglk[ACvsBD] > loglk[ABvsCD] - Constants::closeLogLkLimit) {
                loglk[ACvsBD] = MLQuartetOptimize(*profiles4[0], *profiles4[2], *profiles4[1], *profiles4[3],
                                                  lenACvsBD, /*pStarTest*/ nullptr, &site_likelihoods[ACvsBD * nPos]);
            }
        } else {
            if (options.mlAccuracy > 1 || loglk[ADvsBC] > loglk[ABvsCD] - Constants::closeLogLkLimit) {
                loglk[ADvsBC] = MLQuartetOptimize(*profiles4[0], *profiles4[3], *profiles4[2], *profiles4[1],
                                                  lenADvsBC, /*pStarTest*/ nullptr, &site_likelihoods[ADvsBC * nPos]);
            }
        }

        int choice;
        if (loglk[ABvsCD] >= loglk[ACvsBD] && loglk[ABvsCD] >= loglk[ADvsBC]) {
            choice = ABvsCD;
        } else if (loglk[ACvsBD] >= loglk[ABvsCD] && loglk[ACvsBD] >= loglk[ADvsBC]) {
            choice = ACvsBD;
        } else {
            choice = ADvsBC;
        }
        /* ignore small changes in likelihood */
        bool badSplit = loglk[choice] > loglk[ABvsCD] + Constants::treeLogLkDelta;

        /* constraint penalties, indexed by nni_t (lower is better) */
        double p[3];
        quartetConstraintPenalties(profiles4, /*OUT*/p);
        bool bBadConstr = p[ABvsCD] > p[ACvsBD] + tolerance || p[ABvsCD] > p[ADvsBC] + tolerance;
        bool violateConstraint = false;

        for (int64_t iC = 0; iC < (int64_t) constraintSeqs.size(); iC++) {
            if (splitViolatesConstraint(profiles4, iC)) {
                violateConstraint = true;
                break;
            }
        }
        splitcount.nSplits++;
        if (violateConstraint) {
            splitcount.nConstraintViolations++;
        }
        if (badSplit) {
            splitcount.nBadSplits++;
        }
        if (badSplit && bBadConstr) {
            splitcount.nBadBoth++;
        }
        if (badSplit) {
            double delta = loglk[choice] - loglk[ABvsCD];
            /* If ABvsCD is favored over the more likely NNI by constraints,
           then this is probably a bad split because of the constraint */
            if (p[choice] > p[ABvsCD] + tolerance) {
                splitcount.dWorstDeltaConstrained = std::max(delta, splitcount.dWorstDeltaConstrained);
            } else {
                splitcount.dWorstDeltaUnconstrained = std::max(delta, splitcount.dWorstDeltaUnconstrained);
            }
        }
        if (options.nBootstrap > 0) {
            support[node] = badSplit ? 0.0 : SHSupport(col, loglk, site_likelihoods);
        }

        /* No longer needed */
        upProfiles[nodeABCD[0]].reset();
        upProfiles[nodeABCD[1]].reset();
        upProfiles[nodeABCD[2]].reset();
    }
    return iNodesDone;
}


AbsNeighbourJoining(void)::initNNIStats(std::vector<NNIStats> &stats) {
    stats.resize(maxnode);
    const int64_t LargeAge = 1000000;
    for (int64_t i = 0; i < maxnode; i++) {
        stats[i].delta = 0;
        stats[i].support = 0;
        if (i == root || i < (int64_t) nSeqs) {
            stats[i].age = LargeAge;
            stats[i].subtreeAge = LargeAge;
        } else {
            stats[i].age = 0;
            stats[i].subtreeAge = 0;
        }
    }
}

/* one-dimensional minimization - as input a lower and an upper limit and a trial
  value for the minimum is needed: xmin < xguess < xmax
  the function and a fractional tolerance has to be specified
  onedimenmin returns the optimal x value and the value of the function
  and its second derivative at this point
  */
AbsNeighbourJoining(template<typename Function, typename Data> double)::
onedimenmin(double xmin, double xguess, double xmax, Function f, Data &data, double ftol, double atol, double &fx,
            double &f2x) {
    double optx, ax, bx, cx, fa, fb, fc;

    /* first attempt to bracketize minimum */
    if (xguess == xmin) {
        ax = xmin;
        bx = 2.0 * xguess;
        cx = 10.0 * xguess;
    } else if (xguess <= 2.0 * xmin) {
        ax = xmin;
        bx = xguess;
        cx = 5.0 * xguess;
    } else {
        ax = 0.5 * xguess;
        bx = xguess;
        cx = 2.0 * xguess;
    }
    if (cx > xmax) {
        cx = xmax;
    }
    if (bx >= cx) {
        bx = 0.5 * (ax + cx);
    }
    if (options.verbose > 4) {
        log << strformat("onedimenmin lo %.4f guess %.4f hi %.4f range %.4f %.4f", ax, bx, cx, xmin, xmax) << std::endl;
    }
    /* ideally this range includes the true minimum, i.e.,
       fb < fa and fb < fc
       if not, we gradually expand the boundaries until it does,
       or we near the boundary of the allowed range and use that
    */
    fa = f(ax, data);
    fb = f(bx, data);
    fc = f(cx, data);
    while (fa < fb && ax > xmin) {
        ax = (ax + xmin) / 2.0;
        if (ax < 2.0 * xmin) {    /* give up on shrinking the region */
            ax = xmin;
        }
        fa = f(ax, data);
    }
    while (fc < fb && cx < xmax) {
        cx = (cx + xmax) / 2.0;
        if (cx > xmax * 0.95) {
            cx = xmax;
        }
        fc = f(cx, data);
    }
    optx = brent(ax, bx, cx, f, data, ftol, atol, fx, f2x, fa, fb, fc);

    if (options.verbose > 4) {
        log << strformat("onedimenmin reaches optimum f(%.4f) = %.4f f2x %.4f", optx, fx, f2x) << std::endl;
    }
    return optx; /* return optimal x */
}

/******************************************************************************/
/* Minimization of a 1-dimensional function by Brent's method (Numerical Recipes)
 * Borrowed from Tree-Puzzle 5.1 util.c under GPL
 * Modified by M.N.P to pass in the accessory data for the optimization function,
 * to use 2x bounds around the starting guess and expand them if necessary,
 * and to use both a fractional and an absolute tolerance
 */

#define ITMAX 100
#define CGOLD 0.3819660
#define TINY 1.0e-20
#define ZEPS 1.0e-10
#define SHFT(a, b, c, d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

AbsNeighbourJoining(template<typename Function, typename Data> double)::
brent(double ax, double bx, double cx, Function f, Data &data, double ftol, double atol, double &foptx, double &f2optx,
      double fax, double fbx, double fcx) {
    double a, b, d = 0, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
    double xw, wv, vx;
    double e = 0.0;

    a = (ax < cx ? ax : cx);
    b = (ax > cx ? ax : cx);
    x = bx;
    fx = fbx;
    if (fax < fcx) {
        w = ax;
        fw = fax;
        v = cx;
        fv = fcx;
    } else {
        w = cx;
        fw = fcx;
        v = ax;
        fv = fax;
    }
    for (int64_t iter = 1; iter <= ITMAX; iter++) {
        xm = 0.5 * (a + b);
        tol1 = ftol * std::fabs(x);
        tol2 = 2.0 * (tol1 + ZEPS);
        if (std::fabs(x - xm) <= (tol2 - 0.5 * (b - a)) || std::fabs(a - b) < atol) {
            foptx = fx;
            xw = x - w;
            wv = w - v;
            vx = v - x;
            f2optx = 2.0 * (fv * xw + fx * wv + fw * vx) / (v * v * xw + x * x * wv + w * w * vx);
            return x;
        }
        if (std::fabs(e) > tol1) {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p;
            q = std::fabs(q);
            etemp = e;
            e = d;
            if (std::fabs(p) >= std::fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                d = CGOLD * (e = (x >= xm ? a - x : b - x));
            } else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = SIGN(tol1, xm - x);
            }
        } else {
            d = CGOLD * (e = (x >= xm ? a - x : b - x));
        }
        u = (std::fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
        fu = f(u, data);
        if (fu <= fx) {
            if (u >= x) a = x; else b = x;
            SHFT(v, w, x, u)
            SHFT(fv, fw, fx, fu)
        } else {
            if (u < x) { a = u; } else { b = u; }
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    }
    foptx = fx;
    xw = x - w;
    wv = w - v;
    vx = v - x;
    f2optx = 2.0 * (fv * xw + fx * wv + fw * vx) / (v * v * xw + x * x * wv + w * w * vx);
    return x;
}

#undef ITMAX
#undef CGOLD
#undef TINY
#undef ZEPS
#undef SHFT
#undef SIGN


/* Numerical code for the gamma distribution is modified from the PhyML 3 code
   (GNU public license) of Stephane Guindon
*/

AbsNeighbourJoining(double)::lnGamma(double alpha) {
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
    double x = alpha, f = 0, z;
    if (x < 7) {
        f = 1;
        z = x - 1;
        while (++z < 7) f *= z;
        x = z;
        f = -(double) std::log(f);
    }
    z = 1 / (x * x);
    return f + (x - 0.5) * (double) std::log(x) - x + .918938533204673
           + (((-.000595238095238 * z + .000793650793651) * z - .002777777777778) * z
              + .083333333333333) / x;
}

AbsNeighbourJoining(double)::incompleteGamma(double x, double alpha, double ln_gamma_alpha) {
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
	   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
    int i;
    double p = alpha, g = ln_gamma_alpha;
    double accurate = 1e-8, overflow = 1e30;
    double factor, gin = 0, rn = 0, a = 0, b = 0, an = 0, dif = 0, term = 0, pn[6];

    if (x == 0) return (0);
    if (x < 0 || p <= 0) return (-1);

    factor = (double) std::exp(p * (double) std::log(x) - x - g);
    if (x > 1 && x >= p) goto l30;
    /* (1) series expansion */
    gin = 1;
    term = 1;
    rn = p;
    l20:
    rn++;
    term *= x / rn;
    gin += term;

    if (term > accurate) goto l20;
    gin *= factor / p;
    goto l50;
    l30:
    /* (2) continued fraction */
    a = 1 - p;
    b = a + x + 1;
    term = 0;
    pn[0] = 1;
    pn[1] = x;
    pn[2] = x + 1;
    pn[3] = x * b;
    gin = pn[2] / pn[3];
    l32:
    a++;
    b += 2;
    term++;
    an = a * term;
    for (i = 0; i < 2; i++) pn[i + 4] = b * pn[i + 2] - an * pn[i];
    if (pn[5] == 0) goto l35;
    rn = pn[4] / pn[5];
    dif = fabs(gin - rn);
    if (dif > accurate) goto l34;
    if (dif <= accurate * rn) goto l42;
    l34:
    gin = rn;
    l35:
    for (i = 0; i < 4; i++) pn[i] = pn[i + 2];
    if (fabs(pn[4]) < overflow) goto l32;
    for (i = 0; i < 4; i++) pn[i] /= overflow;
    goto l32;
    l42:
    gin = 1 - factor * gin;

    l50:
    return (gin);
}

AbsNeighbourJoining()::CompareSeeds::CompareSeeds(const std::vector<numeric_t, typename op_t::Allocator> &outDistances,
                                                  const std::vector<int64_t> &compareSeedGaps) :
        outDistances(outDistances), compareSeedGaps(compareSeedGaps) {}

AbsNeighbourJoining(bool)::CompareSeeds::operator()(int64_t seed1, int64_t seed2) const {
    int64_t gapdiff = compareSeedGaps[seed1] - compareSeedGaps[seed2];
    if (gapdiff != 0) {
        return gapdiff < 0; /* fewer gaps is better */
    }
    double outdiff = outDistances[seed1] - outDistances[seed2];
    if (outdiff < 0) {
        return true; /* closer to more nodes is better */
    }
    if (outdiff > 0) {
        return false;
    }
    return true;

}

AbsNeighbourJoining(bool)::CompareHitsByCriterion::operator()(const Besthit &hit1, const Besthit &hit2) const {
    if (hit1.criterion > hit2.criterion) {
        return false;
    }
    return true;
}


AbsNeighbourJoining(bool)::CompareHitsByIJ::operator()(const Besthit &hit1, const Besthit &hit2) const {
    return (hit1.i != hit2.i ? hit1.i - hit2.i : hit1.j - hit2.j) <= 0;
}

#endif