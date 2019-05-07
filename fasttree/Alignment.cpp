
#include "Alignment.h"
#include "Utils.h"
#include <iostream>
#include <sstream>
#include "HashTable.h"

using namespace fasttree;

Alignment::Alignment(const Options &options, std::istream &fp, std::ostream &log) : options(options),
                                                                                    fp(fp),
                                                                                    log(log) {}

void Alignment::readAlignment() {
    std::string buf;

    readline(fp, buf);
    if (buf[0] == '>') {
        /* FASTA, truncate names at any of these */
        auto nameStop = (options.bQuote ? "'\t" : "(),: \t"); /* bQuote supports the -quote option */
        auto seqSkip = " \t";    /* skip these characters in the sequence */
        names.reserve(100);
        seqs.reserve(100);
        nPos = 0;

        do {
            /* loop over lines */
            if (buf[0] == '>') {
                if (!seqs.empty()) {
                    seqs.back().shrink_to_fit();
                }
                /* truncate the name */
                auto found = buf.find_first_of(nameStop);
                if (found != std::string::npos) {
                    buf.resize(found);
                }
                names.push_back(buf);
                names.back().shrink_to_fit();
                seqs.resize(seqs.size() + 1);
            } else {
                /* count non-space characters and append to sequence */
                auto nKeep = buf.find_first_of(seqSkip);
                if (nKeep == std::string::npos) {
                    nKeep = buf.size();
                }
                auto &seq = seqs[names.size() - 1];
                seq.append(buf.begin(), buf.begin() + nKeep);
                if ((int64_t)seq.size() > nPos) {
                    nPos = seq.size();
                }
            }
        } while (readline(fp, buf));

        if (names.size() != seqs.size()) {
            throw std::invalid_argument("No sequence data for last entry " + names.back());
        }

    } else {
        /* PHYLIP interleaved-like format
           Allow arbitrary length names, require spaces between names and sequences
           Allow multiple alignments, either separated by a single empty line (e.g. seqboot output)
           or not.
         */
        int64_t nSeq = 0;
        if (buf.empty()) {
            if (readline(fp, buf)) {
                throw std::invalid_argument("Empty header line followed by EOF");
            }
        }
        std::stringstream aux(buf);
        aux >> nSeq;
        aux >> nPos;
        if (nSeq < 1 || nPos < 1) {
            throw std::invalid_argument("Error parsing header line: " + buf);
        }

        names.resize(nSeq);
        seqs.resize(nSeq);

        int64_t iSeq = 0;
        while (readline(fp, buf)) {
            if (buf.empty() && (iSeq == nSeq || iSeq == 0)) {
                iSeq = 0;
            } else {
                std::string::size_type j = 0; /* character just past end of name */
                if (buf[0] == ' ') {
                    if (names[iSeq].empty()) {
                        throw std::invalid_argument("No name in phylip line: " + buf);
                    }
                } else {
                    j = buf.find_first_of(' ');
                    if (j == std::string::npos || j == 0) {
                        throw std::invalid_argument("No sequence in phylip line: " + buf);
                    }
                    if (iSeq >= nSeq) {
                        throw std::invalid_argument(
                                "No empty line between sequence blocks (is the sequence count wrong?)");
                    }
                    if (names[iSeq].empty()) {
                        names[iSeq] = buf.substr(0, j);
                    } else {
                        /* check the name */
                        if (buf.compare(0, j, names[iSeq]) == 0) {
                            throw std::invalid_argument(
                                    "Wrong name in phylip line " + buf + "\\nExpected " + names[iSeq]);
                        }
                    }
                }
                for (; j < buf.size(); j++) {
                    if (buf[j] != ' ') {
                        if ((int64_t)seqs[iSeq].size() >= nPos) {
                            throw std::invalid_argument(strformat(
                                    "Too many characters (expected %d) for sequence named %s\\nSo far have:\\n%s",
                                    nPos, names[iSeq].c_str(), seqs[iSeq].c_str()));
                        }
                        if (seqs[iSeq].empty()) {
                            seqs[iSeq].reserve(nPos);
                        }
                        seqs[iSeq] += static_cast<char>(std::toupper(buf[j]));
                    }
                }
                if (options.verbose > 10) {
                    log << strformat("Read iSeq %d name %s seqsofar %s",
                                     iSeq, names[iSeq].c_str(), seqs[iSeq].c_str());
                }
                iSeq++;
                if (iSeq == nSeq && (int64_t)seqs[0].size() == nPos) {
                    break; /* finished alignment */
                }
            }/* end else non-empty phylip line */
        }
        if (iSeq != nSeq && iSeq != 0) {
            throw std::invalid_argument(strformat("Wrong number of sequences: expected %d", nSeq));
        }
    }
    /* Check lengths of sequences */
    for (int64_t i = 0; i < (int64_t)seqs.size(); i++) {
        if ((int64_t)seqs[i].size() != nPos) {
            throw std::invalid_argument(strformat(
                    "Wrong number of characters for %s: expected %d but have %d instead.\n"
                    "This sequence may be truncated, or another sequence may be too long.",
                    names[i].c_str(), nPos, seqs[i].size()));
        }
    }
    /* Replace "." with "-" and warn if we find any */
    /* If nucleotide sequences, replace U with T and N with X */
    bool findDot = false;
    #pragma omp parallel for schedule(static), reduction(||:findDot)
    for (int64_t i = 0; i < (int64_t)seqs.size(); i++) {
        for (int64_t j = 0; j < (int64_t)seqs[i].size(); j++) {
            if (seqs[i][j] == '.') {
                seqs[i][j] = '-';
                findDot = true;
            } else if (options.nCodes == 4) {
                if (seqs[i][j] == 'U') {
                    seqs[i][j] = 'T';
                } else if (seqs[i][j] == 'N') {
                    seqs[i][j] = 'X';
                }
            }
        }
    }
    if (findDot) {
        log << "Warning! Found \".\" character(s). These are treated as gaps" << std::endl;
    }

    if (fp.fail() && !fp.eof()) {
        throw std::invalid_argument("Error reading input file");
    }
}

void Alignment::clearAlignmentSeqs() {
    seqs.clear();
}

void Alignment::clearAlignment() {
    names.clear();
    clearAlignmentSeqs();
    nPos = 0;
}

Uniquify::Uniquify(const Alignment &aln) {
    int64_t nUniqueSeq = 0;
    HashTable hashseqs(aln.seqs);

    uniqueSeq.resize(hashseqs.size());     /* iUnique -> seq */
    uniqueFirst.resize(aln.seqs.size());   /* iUnique -> iFirst in aln */
    alnNext.resize(aln.seqs.size(), -1);   /* i in aln -> next, or -1 */
    alnToUniq.resize(aln.seqs.size(), -1); /* i in aln -> iUnique; many -> -1 */

    for (int64_t i = 0; i < (int64_t)aln.seqs.size(); i++) {
        assert(hashseqs.find(aln.seqs[i]) != nullptr);
        int64_t first = hashseqs[aln.seqs[i]];
        if (first == i) {
            /* table uses a pointer to the last repeated sequence as key so if it is affected when move is used,
             * it will not be accessed again */
            uniqueSeq[nUniqueSeq] = std::move(aln.seqs[i]);
            uniqueFirst[nUniqueSeq] = i;
            alnToUniq[i] = nUniqueSeq;
            nUniqueSeq++;
        } else {
            int last = first;
            while (alnNext[last] != -1) {
                last = alnNext[last];
            }
            assert(last >= 0);
            alnNext[last] = i;
            assert(alnToUniq[last] >= 0 && alnToUniq[last] < nUniqueSeq);
            alnToUniq[i] = alnToUniq[last];
        }
    }
    assert((int64_t)nUniqueSeq == (int64_t)uniqueSeq.size());
}
