
#ifndef VERYFASTTREE_DISTANCEMATRIX_H
#define VERYFASTTREE_DISTANCEMATRIX_H

#include "Options.h"

#define MAXCODES 20
namespace veryfasttree {
    template<typename Precision, int Aligment>
    class DistanceMatrix {
    public:
        typedef Precision numeric_t;

        /* Distances between amino acids */
        alignas(Aligment) numeric_t distances[MAXCODES][MAXCODES];

        /* Inverse of the eigenvalue matrix, for rotating a frequency vector
           into eigenspace so that profile similarity computations are
           O(alphabet) not O(alphabet*alphabet) time.
        */
        alignas(Aligment) numeric_t eigeninv[MAXCODES][MAXCODES];
        alignas(Aligment) numeric_t eigenval[MAXCODES];    /* eigenvalues */


        /* eigentot=eigeninv times the all-1s frequency vector
           useful for normalizing rotated frequency vectors
        */
        alignas(Aligment) numeric_t eigentot[MAXCODES];

        /* codeFreq is the transpose of the eigeninv matrix is
           the rotated frequency vector for each code */
        alignas(Aligment) numeric_t codeFreq[MAXCODES][alignsz(MAXCODES, Aligment/ sizeof(numeric_t))];
        alignas(Aligment) numeric_t gapFreq[MAXCODES];

        bool setted;

        void readDistanceMatrix(const Options &options, std::ostream &log);

        void matrixBLOSUM45();

        void setupDistanceMatrix(const Options &options, std::ostream &log);

        operator bool();

    private:

        const static DistanceMatrix<Precision, Aligment> _matrixBLOSUM45;

        void
        readMatrix(const Options &options, const std::string &filename, numeric_t codes[][MAXCODES], bool checkCodes);

        void readVector(const Options &options, const std::string &filename, numeric_t codes[]);

    };

}

#include "DistanceMatrix.tcc"

#endif
