
#ifndef FASTTREE_DISTANCEMATRIX_H
#define FASTTREE_DISTANCEMATRIX_H

#include "Options.h"

#define MAXCODES 20
namespace fasttree {
    template<typename Precision>
    class DistanceMatrix {
    public:
        typedef Precision numeric_t;

        /* Distances between amino acids */
        numeric_t distances[MAXCODES][MAXCODES];

        /* Inverse of the eigenvalue matrix, for rotating a frequency vector
           into eigenspace so that profile similarity computations are
           O(alphabet) not O(alphabet*alphabet) time.
        */
        numeric_t eigeninv[MAXCODES][MAXCODES];
        numeric_t eigenval[MAXCODES];    /* eigenvalues */


        /* eigentot=eigeninv times the all-1s frequency vector
           useful for normalizing rotated frequency vectors
        */
        numeric_t eigentot[MAXCODES];

        /* codeFreq is the transpose of the eigeninv matrix is
           the rotated frequency vector for each code */
        numeric_t codeFreq[MAXCODES][MAXCODES];
        numeric_t gapFreq[MAXCODES];

        bool setted;

        void readDistanceMatrix(const Options &options, std::ostream &log);

        void matrixBLOSUM45();

        void setupDistanceMatrix(const Options &options, std::ostream &log);

        operator bool();

    private:

        const static DistanceMatrix<Precision> _matrixBLOSUM45;

        void readMatrix(const Options &options, const std::string &filename, numeric_t **codes, bool checkCodes);

        void readVector(const Options &options, const std::string &filename, numeric_t *codes);

    };

}

#include "DistanceMatrix.tcc"

#endif
