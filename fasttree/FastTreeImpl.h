
#ifndef FASTTREE_FASTTREEIMPL_H
#define FASTTREE_FASTTREEIMPL_H

#include "Options.h"
#include <memory>
#include <iostream>
#include <fstream>

#define MAXCODES 20

namespace fasttree {
    template<typename Precision, template<typename> typename Operations>
    class FastTreeImpl {
    private:

        typedef Precision numeric_t;

        struct DistanceMatrix {
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
        };

        Options &options;
        std::istream &input;
        std::ostream &output;
        std::ostream &log;

        std::ifstream fpIn;
        std::ifstream fpConstraints;
        std::ifstream fpInTree;

        DistanceMatrix *distanceMatrix;
        DistanceMatrix matrixCustom;
        static DistanceMatrix matrixBLOSUM45;


        void readDistanceMatrix();

        void readMatrix(const std::string &filename, numeric_t codes[MAXCODES][MAXCODES], bool checkCodes);

        void readVector(const std::string &filename, numeric_t codes[MAXCODES]);

        void setupDistanceMatrix();

        void readAlignment(std::istream &fp, bool bQuote);


    public:
        FastTreeImpl(Options &options, std::istream &input, std::ostream &output, std::ostream &log);

        void run();
    };
}

#include "FastTreeImpl.tcc"

#endif
