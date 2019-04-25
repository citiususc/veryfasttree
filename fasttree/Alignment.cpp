
#include "Alignment.h"
#include "Utils.h"
#include <iostream>
#include <sstream>

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
                _nSeq++;
            } else {
                /* count non-space characters and append to sequence */
                auto nKeep = buf.find_first_of(seqSkip);
                if (nKeep == std::string::npos) {
                    nKeep = buf.size();
                }
                auto &seq = seqs[names.size() - 1];
                seq.append(buf.begin(), buf.begin() + nKeep);
                if (seq.size() > _nPos) {
                    _nPos = seq.size();
                }
            }
        } while (readline(fp, buf));

        if (names.size() != seqs.size()) {
            throw new std::invalid_argument("No sequence data for last entry " + names.back());
        }

    } else {
        /* PHYLIP interleaved-like format
           Allow arbitrary length names, require spaces between names and sequences
           Allow multiple alignments, either separated by a single empty line (e.g. seqboot output)
           or not.
         */
        if (buf.empty()) {
            if (readline(fp, buf)) {
                throw new std::invalid_argument("Empty header line followed by EOF");
            }
        }
        std::stringstream aux(buf);
        aux >> _nSeq;
        aux >> _nPos;
        if (_nSeq < 1 || _nPos < 1) {
            throw new std::invalid_argument("Error parsing header line: " + buf);
        }

        for (size_t i = 0; i < _nSeq; i++) {
            names[i] = "";
            seqs[i].reserve(_nPos);
        }

        size_t iSeq = 0;
        while (readline(fp, buf)) {
            if (buf.empty() && (iSeq == _nSeq || iSeq == 0)) {
                iSeq = 0;
            } else {
                size_t j = 0; /* character just past end of name */
                if (buf[0] == ' ') {
                    if (names[iSeq].empty()) {
                        throw new std::invalid_argument("No name in phylip line: " + buf);
                    }
                } else {
                    j = buf.find_first_of(' ');
                    if (j != std::string::npos || j == 0) {
                        throw new std::invalid_argument("No sequence in phylip line: " + buf);
                    }
                    if (iSeq >= _nSeq) {
                        throw new std::invalid_argument(
                                "No empty line between sequence blocks (is the sequence count wrong?)");
                    }
                    if (names[iSeq].empty()) {
                        names[iSeq] = buf.substr(0, j);
                    } else {
                        /* check the name */
                        if (buf.compare(0, j, names[iSeq]) == 0) {
                            throw new std::invalid_argument(
                                    "Wrong name in phylip line " + buf + "\\nExpected " + names[iSeq]);
                        }
                    }
                }
                for (; j < buf.size(); j++) {
                    if (buf[j] != ' ') {
                        if (seqs[iSeq].size() >= _nPos) {
                            throw new std::invalid_argument(strformat(
                                    "Too many characters (expected %d) for sequence named %s\\nSo far have:\\n%s",
                                    _nPos, names[iSeq].c_str(), seqs[iSeq].c_str()));
                        }
                        seqs[iSeq] += std::toupper(buf[j]);
                    }
                }
                if (options.verbose > 10) {
                    log << strformat("Read iSeq %d name %s seqsofar %s",
                                     iSeq, names[iSeq].c_str(), seqs[iSeq].c_str());
                }
                iSeq++;
                if (iSeq == _nSeq && seqs[0].size() == _nPos) {
                    break; /* finished alignment */
                }
            }/* end else non-empty phylip line */
        }
        if (iSeq != _nSeq && iSeq != 0) {
            throw new std::invalid_argument(strformat("Wrong number of sequences: expected %d", _nSeq));
        }
    }
    /* Check lengths of sequences */
    for (size_t i = 0; i < _nSeq; i++) {
        if (seqs[i].size() != _nPos) {
            throw new std::invalid_argument(strformat(
                    "Wrong number of characters for %s: expected %d but have %d instead.\n"
                    "This sequence may be truncated, or another sequence may be too long.",
                    names[i], _nPos, seqs[i].size()));
        }
    }
    /* Replace "." with "-" and warn if we find any */
    /* If nucleotide sequences, replace U with T and N with X */
    bool findDot = false;
    for (size_t i = 0; i < _nSeq; i++) {//TODO OpenMP
        for (size_t j = 0; j < seqs[i].size(); j++) {
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
        throw new std::invalid_argument("Error reading input file");
    }
}

void Alignment::clearAlignmentSeqs() {
    seqs.clear();
}

void Alignment::clearAlignment() {
    names.clear();
    clearAlignmentSeqs();
    _nSeq = 0;
    _nPos = 0;
}

const size_t &Alignment::nPos() const {
    return _nPos;
}

const size_t &Alignment::nSeq() const {
    return _nSeq;
}
