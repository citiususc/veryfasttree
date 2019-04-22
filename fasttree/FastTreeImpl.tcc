
#ifndef FASTTREE_FASTTREEIMPL_TCC
#define FASTTREE_FASTTREEIMPL_TCC

#include "FastTreeImpl.h"
#include "Utils.h"

#define AbsFastTreeImpl \
template<typename Precision, template<typename> typename Operations, template<typename> typename Allocator> \
FastTreeImpl<Precision, Operations, Allocator>


namespace fasttree {

    AbsFastTreeImpl::FastTreeImpl(Options &options) : options(options) {


        /**
        if (options.matrixPrefix.size()>0) {
            if (!options.useMatrix) {
                throw std::invalid_argument("Cannot use both -matrix and -nomatrix arguments!");
            }
            distance_matrix = ReadDistanceMatrix(matrixPrefix);
        } else if (options.useMatrix) {   // use default matrix
            assert(options.nCodes == 20);
            distance_matrix = &matrixBLOSUM45;
            SetupDistanceMatrix(distance_matrix);
        } else {
            distance_matrix = NULL;
        }

        int iAln;
        FILE *fpIn = fileName != NULL ? fopen(fileName, "r") : stdin;
        if (fpIn == NULL) {
            fprintf(stderr, "Cannot read %s\n", fileName);
            exit(1);
        }
        FILE *fpConstraints = NULL;
        if (constraintsFile != NULL) {
            fpConstraints = fopen(constraintsFile, "r");
            if (fpConstraints == NULL) {
                fprintf(stderr, "Cannot read %s\n", constraintsFile);
                exit(1);
            }
        }

        FILE *fpInTree = NULL;
        if (intreeFile != NULL) {
            fpInTree = fopen(intreeFile, "r");
            if (fpInTree == NULL) {
                fprintf(stderr, "Cannot read %s\n", intreeFile);
                exit(1);
            }
        }**/

    }
}
#endif