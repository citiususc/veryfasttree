
#include "TransitionMatrix.h"

#define AbsTransitionMatrix(...) \
template<typename Precision, int Aligment> \
__VA_ARGS__ veryfasttree::TransitionMatrix<Precision, Aligment>

#define double_alignas alignas(Aligment < sizeof(double) ? sizeof(double) : Aligment)

AbsTransitionMatrix()::TransitionMatrix() : setted(false) {}

AbsTransitionMatrix()::operator bool() { return setted; }

AbsTransitionMatrix(void)::createTransitionMatrixJTT92(const Options &options) {
    createTransitionMatrix(options, matrixJTT92, statJTT92);
}

AbsTransitionMatrix(void)::createTransitionMatrixWAG01(const Options &options) {
    createTransitionMatrix(options, matrixWAG01, statWAG01);
}

AbsTransitionMatrix(void)::createTransitionMatrixLG08(const Options &options) {
    createTransitionMatrix(options, matrixLG08, statLG08);
}

AbsTransitionMatrix(void)::createGTR(const Options &options, double r[], double f[]) {
    double matrix[4][MAXCODES];
    assert(options.nCodes == 4);


    /* Place rates onto a symmetric matrix, but correct by f(target), so that
       stationary distribution f[] is maintained
       Leave diagonals as 0 (CreateTransitionMatrix will fix them)
    */
    int imat = 0;
    for (int i = 0; i < options.nCodes; i++) {
        matrix[i][i] = 0;
        for (int j = i + 1; j < options.nCodes; j++) {
            double rate = r[imat++];
            assert(rate > 0);
            /* Want t(matrix) * f to be 0 */
            matrix[i][j] = rate * f[i];
            matrix[j][i] = rate * f[j];
        }
    }
    /* Compute average mutation rate */
    double total_rate = 0;
    for (int i = 0; i < options.nCodes; i++) {
        for (int j = 0; j < options.nCodes; j++) {
            total_rate += f[i] * matrix[i][j];
        }
    }
    assert(total_rate > 1e-6);
    double inv = 1.0 / total_rate;
    for (int i = 0; i < options.nCodes; i++) {
        for (int j = 0; j < options.nCodes; j++) {
            matrix[i][j] *= inv;
        }
    }
    createTransitionMatrix(options, matrix, f);
}

AbsTransitionMatrix(void)::readAATransitionMatrix(const Options &options, /*IN*/ const std::string &filename) {
    assert(options.nCodes == 20);
    double stat[20];
    double matrix[MAXCODES][MAXCODES];
    std::string buf;
    std::ifstream fp(filename);
    if (fp.fail()) {
        throw std::invalid_argument("Cannot read transition matrix file " + filename);
    }
    std::string expected;
    expected.reserve(2 * MAXCODES + 20);
    for (int i = 0; i < 20; i++) {
        expected.insert(expected.end(), Constants::codesStringAA[i]);
        expected.insert(expected.end(), '\t');
    }
    expected.insert(expected.end(), '*');

    if (!readline(fp, buf)) {
        throw std::invalid_argument("Error reading header line from transition matrix file");
    }
    if (buf != expected) {
        throw std::invalid_argument("Invalid header line in transition matrix file, it must match: " + expected);
    }
    for (int i = 0; i < 20; i++) {
        if (!readline(fp, buf)) {
            throw std::invalid_argument("Error reading matrix line");
        }
        std::istringstream sfields(buf);
        std::string field;

        if (!std::getline(sfields, field, '\t') || field.size() != 1 || field[0] != Constants::codesStringAA[i]) {
            throw std::invalid_argument(strformat("Line for amino acid %c does not have the expected beginning",
                                                  Constants::codesStringAA[i]));
        }
        for (int j = 0; j < 20; j++) {
            if (!std::getline(sfields, field, '\t')) {
                throw std::invalid_argument(
                        strformat("Not enough fields for amino acid %c", Constants::codesStringAA[i]));
            }
            matrix[i][j] = std::stod(field);
        }
        if (!std::getline(sfields, field, '\t')) {
            throw std::invalid_argument(strformat("Not enough fields for amino acid %c", Constants::codesStringAA[i]));
        }
        stat[i] = std::stod(field);
    }

    double tol = 1e-5;
    /* Verify that stat is positive and sums to 1 */
    double statTot = 0;
    for (int i = 0; i < 20; i++) {
        if (stat[i] < tol) {
            throw std::invalid_argument(strformat("stationary frequency for amino acid %c must be positive",
                                                  Constants::codesStringAA[i]));
        }
        statTot += stat[i];
    }
    if (fabs(statTot - 1) > tol) {
        throw std::invalid_argument(strformat("stationary frequencies must sum to 1 -- actual sum is %g", statTot));
    }

    /* Verify that diagonals are negative and dot product of stat and diagonals is -1 */
    double totRate = 0;
    for (int i = 0; i < 20; i++) {
        double diag = matrix[i][i];
        if (diag > -tol) {
            throw std::invalid_argument(strformat("transition rate(%c,%c) must be negative",
                                                  Constants::codesStringAA[i], Constants::codesStringAA[i]));
        }
        totRate += stat[i] * diag;
    }
    if (fabs(totRate + 1) > tol) {
        throw std::invalid_argument(strformat("Dot product of matrix diagonal and stationary frequencies must "
                                              "be -1 -- actual dot product is %g", totRate));
    }

    /* Verify that each off-diagonal entry is nonnegative and that each column sums to 0 */
    for (int j = 0; j < 20; j++) {
        double colSum = 0;
        for (int i = 0; i < 20; i++) {
            double value = matrix[i][j];
            colSum += value;
            if (i != j && value < 0) {
                throw std::invalid_argument(strformat("Off-diagonal matrix entry for (%c,%c) is negative",
                                                      Constants::codesStringAA[i], Constants::codesStringAA[j]));
            }
        }
        if (fabs(colSum) > tol) {
            throw std::invalid_argument(strformat("Sum of column %c must be zero -- actual sum is %g",
                                                  Constants::codesStringAA[j], colSum));
        }
    }
    return createTransitionMatrix(options, matrix, stat);
}

AbsTransitionMatrix(void)::
createTransitionMatrix(const Options &options, const double matrix[MAXCODES][MAXCODES], const double _stat[MAXCODES]) {
    double_alignas double sqrtstat[20];
    auto nCodes = options.nCodes;

    for (int i = 0; i < options.nCodes; i++) {
        stat[i] = _stat[i];
        statinv[i] = 1.0 / _stat[i];
        sqrtstat[i] = std::sqrt(_stat[i]);
    }

    double_alignas double sym[20 * 20];        /* symmetrized matrix M' */
    /* set diagonals so columns sums are 0 before symmetrization */
    for (int i = 0; i < nCodes; i++) {
        for (int j = 0; j < nCodes; j++) {
            sym[nCodes * i + j] = matrix[i][j];
        }
    }
    for (int j = 0; j < nCodes; j++) {
        double sum = 0;
        sym[nCodes * j + j] = 0;
        for (int i = 0; i < nCodes; i++) {
            sum += sym[nCodes * i + j];
        }
        sym[nCodes * j + j] = -sum;
    }
    /* M' = S**-1 M S */
    for (int i = 0; i < nCodes; i++) {
        for (int j = 0; j < nCodes; j++) {
            sym[nCodes * i + j] *= sqrtstat[j] / sqrtstat[i];
        }
    }

    /* eigen decomposition of M' -- note that eigenW is the transpose of what we want,
       which is eigenvectors in columns */
    double eigenW[20 * 20], eval[20], e[20];
    for (int i = 0; i < nCodes * nCodes; i++) {
        eigenW[i] = sym[i];
    }
    tred2(eigenW, nCodes, nCodes, eval, e);
    tqli(eval, e, nCodes, nCodes, eigenW);

    /* save eigenvalues */
    for (int i = 0; i < nCodes; i++) {
        this->eigenval[i] = eval[i];
    }

    /* compute eigen decomposition of M into t(codeFreq): V = S*W */
    /* compute inverse of V in eigeninv: V**-1 = t(W) S**-1  */
    for (int i = 0; i < nCodes; i++) {
        for (int j = 0; j < nCodes; j++) {
            this->eigeninv[i][j] = eigenW[nCodes * i + j] / sqrtstat[j];
            this->eigeninvT[j][i] = this->eigeninv[i][j];
        }
    }
    for (int i = 0; i < nCodes; i++) {
        for (int j = 0; j < nCodes; j++) {
            this->codeFreq[i][j] = eigenW[j * nCodes + i] * sqrtstat[i];
        }
    }
    /* codeFreq[NOCODE] is the rotation of (1,1,...) not (1/nCodes,1/nCodes,...), which
       gives correct posterior probabilities
    */
    for (int j = 0; j < nCodes; j++) {
        this->codeFreq[NOCODE][j] = 0.0;
        for (int i = 0; i < nCodes; i++) {
            this->codeFreq[NOCODE][j] += this->codeFreq[i][j];
        }
    }
    /* save some posterior probabilities for approximating later:
       first, we compute P(B | A, t) for t = approxMLnearT, by using
       V * exp(L*t) * V**-1 */
    double_alignas double expvalues[MAXCODES];
    for (int i = 0; i < nCodes; i++) {
        expvalues[i] = std::exp(Constants::approxMLnearT * this->eigenval[i]);
    }
    double_alignas double LVinv[MAXCODES][MAXCODES]; /* exp(L*t) * V**-1 */
    for (int i = 0; i < nCodes; i++) {
        for (int j = 0; j < nCodes; j++) {
            LVinv[i][j] = this->eigeninv[i][j] * expvalues[i];
        }
    }
    /* matrix transform for converting A -> B given t: transt[i][j] = P(j->i | t) */
    double_alignas double transt[MAXCODES][MAXCODES];
    for (int i = 0; i < nCodes; i++) {
        for (int j = 0; j < nCodes; j++) {
            transt[i][j] = 0;
            for (int k = 0; k < nCodes; k++) {
                transt[i][j] += this->codeFreq[i][k] * LVinv[k][j];
            }
        }
    }
    /* nearP[i][j] = P(parent = j | both children are i) = P(j | i,i) ~ stat(j) * P(j->i | t)**2 */
    for (int i = 0; i < nCodes; i++) {
        double_alignas double _nearP[MAXCODES];
        double tot = 0;
        for (int j = 0; j < nCodes; j++) {
            assert(transt[j][i] > 0);
            assert(this->stat[j] > 0);
            _nearP[j] = this->stat[j] * transt[i][j] * transt[i][j];
            tot += _nearP[j];
        }
        assert(tot > 0);
        for (int j = 0; j < nCodes; j++) {
            _nearP[j] *= 1.0 / tot;
        }
        /* save nearP in transmat->nearP[i][] */
        for (int j = 0; j < nCodes; j++) {
            this->nearP[i][j] = _nearP[j];
        }
        /* multiply by 1/stat and rotate nearP */
        for (int j = 0; j < nCodes; j++) {
            _nearP[j] /= this->stat[j];
        }
        for (int j = 0; j < nCodes; j++) {
            double rot = 0;
            for (int k = 0; k < nCodes; k++) {
                rot += _nearP[k] * this->codeFreq[i][j];
            }
            this->nearFreq[i][j] = rot;
        }
    }
    setted = true;
}

// @formatter:off
AbsTransitionMatrix(const double)::statJTT92[MAXCODES] = {0.07674789,0.05169087,0.04264509,0.05154407,0.01980301,0.04075195,0.06182989,0.07315199,0.02294399,0.05376110,0.09190390,0.05867583,0.02382594,0.04012589,0.05090097,0.06876503,0.05856501,0.01426057,0.03210196,0.06600504};
AbsTransitionMatrix(const double)::matrixJTT92[MAXCODES][MAXCODES] = {
        { -1.247831,0.044229,0.041179,0.061769,0.042704,0.043467,0.08007,0.136501,0.02059,0.027453,0.022877,0.02669,0.041179,0.011439,0.14794,0.288253,0.362223,0.006863,0.008388,0.227247 },
        { 0.029789,-1.025965,0.023112,0.008218,0.058038,0.159218,0.014895,0.070364,0.168463,0.011299,0.019517,0.33179,0.022599,0.002568,0.038007,0.051874,0.032871,0.064714,0.010272,0.008731 },
        { 0.022881,0.019068,-1.280568,0.223727,0.014407,0.03644,0.024576,0.034322,0.165676,0.019915,0.005085,0.11144,0.012712,0.004237,0.006356,0.213134,0.098304,0.00339,0.029661,0.00678 },
        { 0.041484,0.008194,0.270413,-1.044903,0.005121,0.025095,0.392816,0.066579,0.05736,0.005634,0.003585,0.013316,0.007682,0.002049,0.007682,0.030217,0.019462,0.002049,0.023559,0.015877 },
        { 0.011019,0.022234,0.00669,0.001968,-0.56571,0.001771,0.000984,0.011609,0.013577,0.003345,0.004526,0.001377,0.0061,0.015348,0.002755,0.043878,0.008264,0.022628,0.041124,0.012199 },
        { 0.02308,0.125524,0.034823,0.019841,0.003644,-1.04415,0.130788,0.010528,0.241735,0.003644,0.029154,0.118235,0.017411,0.00162,0.066406,0.021461,0.020651,0.007288,0.009718,0.008098 },
        { 0.064507,0.017816,0.035632,0.471205,0.003072,0.198435,-0.944343,0.073107,0.015973,0.007372,0.005529,0.111197,0.011058,0.003072,0.011058,0.01843,0.019659,0.006143,0.0043,0.027646 },
        { 0.130105,0.099578,0.058874,0.09449,0.042884,0.018898,0.086495,-0.647831,0.016717,0.004361,0.004361,0.019625,0.010176,0.003634,0.017444,0.146096,0.023986,0.039976,0.005815,0.034162 },
        { 0.006155,0.074775,0.089138,0.025533,0.01573,0.1361,0.005927,0.005243,-1.135695,0.003648,0.012767,0.010259,0.007523,0.009119,0.026217,0.016642,0.010487,0.001824,0.130629,0.002508 },
        { 0.01923,0.011752,0.025106,0.005876,0.009081,0.004808,0.00641,0.003205,0.008547,-1.273602,0.122326,0.011218,0.25587,0.047542,0.005342,0.021367,0.130873,0.004808,0.017094,0.513342 },
        { 0.027395,0.0347,0.010958,0.006392,0.021003,0.065748,0.008219,0.005479,0.051137,0.209115,-0.668139,0.012784,0.354309,0.226465,0.093143,0.053877,0.022829,0.047485,0.021916,0.16437 },
        { 0.020405,0.376625,0.153332,0.015158,0.004081,0.170239,0.105525,0.015741,0.026235,0.012243,0.008162,-0.900734,0.037896,0.002332,0.012243,0.027401,0.06005,0.00583,0.004664,0.008162 },
        { 0.012784,0.010416,0.007102,0.003551,0.007339,0.01018,0.004261,0.003314,0.007812,0.113397,0.091854,0.015388,-1.182051,0.01018,0.003788,0.006865,0.053503,0.005682,0.004261,0.076466 },
        { 0.00598,0.001993,0.003987,0.001595,0.031098,0.001595,0.001993,0.001993,0.015948,0.035484,0.098877,0.001595,0.017144,-0.637182,0.006778,0.03668,0.004784,0.021131,0.213701,0.024719 },
        { 0.098117,0.037426,0.007586,0.007586,0.007081,0.082944,0.009104,0.012138,0.058162,0.005058,0.051587,0.010621,0.008092,0.008598,-0.727675,0.144141,0.059679,0.003035,0.005058,0.011632 },
        { 0.258271,0.069009,0.343678,0.040312,0.152366,0.036213,0.020498,0.137334,0.049878,0.02733,0.040312,0.032113,0.019814,0.06286,0.194728,-1.447863,0.325913,0.023914,0.043045,0.025964 },
        { 0.276406,0.037242,0.135003,0.022112,0.02444,0.029677,0.018621,0.019203,0.026768,0.142567,0.014548,0.059936,0.131511,0.006983,0.068665,0.27757,-1.335389,0.006983,0.01222,0.065174 },
        { 0.001275,0.017854,0.001134,0.000567,0.016295,0.002551,0.001417,0.007793,0.001134,0.001275,0.007368,0.001417,0.003401,0.00751,0.00085,0.004959,0.0017,-0.312785,0.010061,0.003542 },
        { 0.003509,0.006379,0.022328,0.014673,0.066664,0.007655,0.002233,0.002552,0.182769,0.010207,0.007655,0.002552,0.005741,0.170967,0.00319,0.020095,0.006698,0.022647,-0.605978,0.005103 },
        { 0.195438,0.011149,0.010493,0.020331,0.040662,0.013117,0.029512,0.030824,0.007214,0.630254,0.11805,0.009182,0.211834,0.040662,0.015084,0.024922,0.073453,0.016396,0.010493,-1.241722 }
};

AbsTransitionMatrix(const double)::statWAG01[MAXCODES] = {0.0866279,0.043972, 0.0390894,0.0570451,0.0193078,0.0367281,0.0580589,0.0832518,0.0244314,0.048466, 0.086209, 0.0620286,0.0195027,0.0384319,0.0457631,0.0695179,0.0610127,0.0143859,0.0352742,0.0708956};
AbsTransitionMatrix(const double)::matrixWAG01[MAXCODES][MAXCODES] = {
        {-1.117151, 0.050147, 0.046354, 0.067188, 0.093376, 0.082607, 0.143908, 0.128804, 0.028817, 0.017577, 0.036177, 0.082395, 0.081234, 0.019138, 0.130789, 0.306463, 0.192846, 0.010286, 0.021887, 0.182381},
        {0.025455, -0.974318, 0.029321, 0.006798, 0.024376, 0.140086, 0.020267, 0.026982, 0.098628, 0.008629, 0.022967, 0.246964, 0.031527, 0.004740, 0.031358, 0.056495, 0.025586, 0.053714, 0.017607, 0.011623},
        {0.020916, 0.026065, -1.452438, 0.222741, 0.010882, 0.063328, 0.038859, 0.046176, 0.162306, 0.022737, 0.005396, 0.123567, 0.008132, 0.003945, 0.008003, 0.163042, 0.083283, 0.002950, 0.044553, 0.008051},
        {0.044244, 0.008819, 0.325058, -0.989665, 0.001814, 0.036927, 0.369645, 0.051822, 0.055719, 0.002361, 0.005077, 0.028729, 0.006212, 0.002798, 0.025384, 0.064166, 0.022443, 0.007769, 0.019500, 0.009120},
        {0.020812, 0.010703, 0.005375, 0.000614, -0.487357, 0.002002, 0.000433, 0.006214, 0.005045, 0.003448, 0.007787, 0.001500, 0.007913, 0.008065, 0.002217, 0.028525, 0.010395, 0.014531, 0.011020, 0.020307},
        {0.035023, 0.117008, 0.059502, 0.023775, 0.003809, -1.379785, 0.210830, 0.012722, 0.165524, 0.004391, 0.033516, 0.150135, 0.059565, 0.003852, 0.035978, 0.039660, 0.033070, 0.008316, 0.008777, 0.011613},
        {0.096449, 0.026759, 0.057716, 0.376214, 0.001301, 0.333275, -1.236894, 0.034593, 0.034734, 0.007763, 0.009400, 0.157479, 0.019202, 0.004944, 0.041578, 0.042955, 0.050134, 0.009540, 0.011961, 0.035874},
        {0.123784, 0.051085, 0.098345, 0.075630, 0.026795, 0.028838, 0.049604, -0.497615, 0.021792, 0.002661, 0.005356, 0.032639, 0.015212, 0.004363, 0.021282, 0.117240, 0.019732, 0.029444, 0.009052, 0.016361},
        {0.008127, 0.054799, 0.101443, 0.023863, 0.006384, 0.110105, 0.014616, 0.006395, -0.992342, 0.003543, 0.012807, 0.022832, 0.010363, 0.017420, 0.017851, 0.018979, 0.012136, 0.006733, 0.099319, 0.003035},
        {0.009834, 0.009511, 0.028192, 0.002006, 0.008654, 0.005794, 0.006480, 0.001549, 0.007029, -1.233162, 0.161294, 0.016472, 0.216559, 0.053891, 0.005083, 0.016249, 0.074170, 0.010808, 0.021372, 0.397837},
        {0.036002, 0.045028, 0.011900, 0.007673, 0.034769, 0.078669, 0.013957, 0.005547, 0.045190, 0.286902, -0.726011, 0.023303, 0.439180, 0.191376, 0.037625, 0.031191, 0.029552, 0.060196, 0.036066, 0.162890},
        {0.058998, 0.348377, 0.196082, 0.031239, 0.004820, 0.253558, 0.168246, 0.024319, 0.057967, 0.021081, 0.016767, -1.124580, 0.060821, 0.005783, 0.036254, 0.062960, 0.090292, 0.008952, 0.008675, 0.019884},
        {0.018288, 0.013983, 0.004057, 0.002124, 0.007993, 0.031629, 0.006450, 0.003564, 0.008272, 0.087143, 0.099354, 0.019123, -1.322098, 0.024370, 0.003507, 0.010109, 0.031033, 0.010556, 0.008769, 0.042133},
        {0.008490, 0.004143, 0.003879, 0.001885, 0.016054, 0.004030, 0.003273, 0.002014, 0.027402, 0.042734, 0.085315, 0.003583, 0.048024, -0.713669, 0.006512, 0.022020, 0.006934, 0.061698, 0.260332, 0.026213},
        {0.069092, 0.032635, 0.009370, 0.020364, 0.005255, 0.044829, 0.032773, 0.011698, 0.033438, 0.004799, 0.019973, 0.026747, 0.008229, 0.007754, -0.605590, 0.077484, 0.038202, 0.006695, 0.010376, 0.015124},
        {0.245933, 0.089317, 0.289960, 0.078196, 0.102703, 0.075066, 0.051432, 0.097899, 0.054003, 0.023306, 0.025152, 0.070562, 0.036035, 0.039831, 0.117705, -1.392239, 0.319421, 0.038212, 0.057419, 0.016981},
        {0.135823, 0.035501, 0.129992, 0.024004, 0.032848, 0.054936, 0.052685, 0.014461, 0.030308, 0.093371, 0.020915, 0.088814, 0.097083, 0.011008, 0.050931, 0.280341, -1.154973, 0.007099, 0.018643, 0.088894},
        {0.001708, 0.017573, 0.001086, 0.001959, 0.010826, 0.003257, 0.002364, 0.005088, 0.003964, 0.003208, 0.010045, 0.002076, 0.007786, 0.023095, 0.002105, 0.007908, 0.001674, -0.466694, 0.037525, 0.005516},
        {0.008912, 0.014125, 0.040205, 0.012058, 0.020133, 0.008430, 0.007267, 0.003836, 0.143398, 0.015555, 0.014757, 0.004934, 0.015861, 0.238943, 0.007998, 0.029135, 0.010779, 0.092011, -0.726275, 0.011652},
        {0.149259, 0.018739, 0.014602, 0.011335, 0.074565, 0.022417, 0.043805, 0.013932, 0.008807, 0.581952, 0.133956, 0.022726, 0.153161, 0.048356, 0.023429, 0.017317, 0.103293, 0.027186, 0.023418, -1.085487},
};

/* Le-Gascuel 2008 model data from Harry Yoo
   https://github.com/hyoo/FastTree
*/
AbsTransitionMatrix(const double)::statLG08[MAXCODES] = {0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767, 0.071586, 0.057337, 0.022355, 0.062157, 0.099081, 0.0646, 0.022951, 0.042302, 0.04404, 0.061197, 0.053287, 0.012066, 0.034155, 0.069147};
AbsTransitionMatrix(const double)::matrixLG08[MAXCODES][MAXCODES] = {
        {-1.08959879,0.03361031,0.02188683,0.03124237,0.19680136,0.07668542,0.08211337,0.16335306,0.02837339,0.01184642,0.03125763,0.04242021,0.08887270,0.02005907,0.09311189,0.37375830,0.16916131,0.01428853,0.01731216,0.20144931},
        {0.02378006,-0.88334349,0.04206069,0.00693409,0.02990323,0.15707674,0.02036079,0.02182767,0.13574610,0.00710398,0.01688563,0.35388551,0.02708281,0.00294931,0.01860218,0.04800569,0.03238902,0.03320688,0.01759004,0.00955956},
        {0.01161996,0.03156149,-1.18705869,0.21308090,0.02219603,0.07118238,0.02273938,0.06034785,0.18928374,0.00803870,0.00287235,0.09004368,0.01557359,0.00375798,0.00679131,0.16825837,0.08398226,0.00190474,0.02569090,0.00351296},
        {0.02096312,0.00657599,0.26929909,-0.86328733,0.00331871,0.02776660,0.27819699,0.04482489,0.04918511,0.00056712,0.00079981,0.01501150,0.00135537,0.00092395,0.02092662,0.06579888,0.02259266,0.00158572,0.00716768,0.00201422},
        {0.03220119,0.00691547,0.00684065,0.00080928,-0.86781864,0.00109716,0.00004527,0.00736456,0.00828668,0.00414794,0.00768465,0.00017162,0.01156150,0.01429859,0.00097521,0.03602269,0.01479316,0.00866942,0.01507844,0.02534728},
        {0.03953956,0.11446966,0.06913053,0.02133682,0.00345736,-1.24953177,0.16830979,0.01092385,0.19623161,0.00297003,0.02374496,0.13185209,0.06818543,0.00146170,0.02545052,0.04989165,0.04403378,0.00962910,0.01049079,0.00857458},
        {0.07434507,0.02605508,0.03877888,0.37538659,0.00025048,0.29554848,-0.84254259,0.02497249,0.03034386,0.00316875,0.00498760,0.12936820,0.01243696,0.00134660,0.03002373,0.04380857,0.04327684,0.00557310,0.00859294,0.01754095},
        {0.11846020,0.02237238,0.08243001,0.04844538,0.03263985,0.01536392,0.02000178,-0.50414422,0.01785951,0.00049912,0.00253779,0.01700817,0.00800067,0.00513658,0.01129312,0.09976552,0.00744439,0.01539442,0.00313512,0.00439779},
        {0.00802225,0.05424651,0.10080372,0.02072557,0.01431930,0.10760560,0.00947583,0.00696321,-1.09324335,0.00243405,0.00818899,0.01558729,0.00989143,0.01524917,0.01137533,0.02213166,0.01306114,0.01334710,0.11863394,0.00266053},
        {0.00931296,0.00789336,0.01190322,0.00066446,0.01992916,0.00452837,0.00275137,0.00054108,0.00676776,-1.41499789,0.25764421,0.00988722,0.26563382,0.06916358,0.00486570,0.00398456,0.06425393,0.00694043,0.01445289,0.66191466},
        {0.03917027,0.02990732,0.00677980,0.00149374,0.05885464,0.05771026,0.00690325,0.00438541,0.03629495,0.41069624,-0.79375308,0.01362360,0.62543296,0.25688578,0.02467704,0.01806113,0.03001512,0.06139358,0.02968934,0.16870919},
        {0.03465896,0.40866276,0.13857164,0.01827910,0.00085698,0.20893479,0.11674330,0.01916263,0.04504313,0.01027583,0.00888247,-0.97644156,0.04241650,0.00154510,0.02521473,0.04836478,0.07344114,0.00322392,0.00852278,0.01196402},
        {0.02579765,0.01111131,0.00851489,0.00058635,0.02051079,0.03838702,0.00398738,0.00320253,0.01015515,0.09808327,0.14487451,0.01506968,-1.54195698,0.04128536,0.00229163,0.00796306,0.04636929,0.01597787,0.01104642,0.04357735},
        {0.01073203,0.00223024,0.00378708,0.00073673,0.04675419,0.00151673,0.00079574,0.00378966,0.02885576,0.04707045,0.10967574,0.00101178,0.07609486,-0.81061579,0.00399600,0.01530562,0.00697985,0.10394083,0.33011973,0.02769432},
        {0.05186360,0.01464471,0.00712508,0.01737179,0.00331981,0.02749383,0.01847072,0.00867414,0.02240973,0.00344749,0.01096857,0.01718973,0.00439734,0.00416018,-0.41664685,0.05893117,0.02516738,0.00418956,0.00394655,0.01305787},
        {0.28928853,0.05251612,0.24529879,0.07590089,0.17040121,0.07489439,0.03745080,0.10648187,0.06058559,0.00392302,0.01115539,0.04581702,0.02123285,0.02214217,0.08188943,-1.42842431,0.39608294,0.01522956,0.02451220,0.00601987},
        {0.11400727,0.03085239,0.10660988,0.02269274,0.06093244,0.05755704,0.03221430,0.00691855,0.03113348,0.05508469,0.01614250,0.06057985,0.10765893,0.00879238,0.03045173,0.34488735,-1.23444419,0.00750412,0.01310009,0.11660005},
        {0.00218053,0.00716244,0.00054751,0.00036065,0.00808574,0.00284997,0.00093936,0.00323960,0.00720403,0.00134729,0.00747646,0.00060216,0.00840002,0.02964754,0.00114785,0.00300276,0.00169919,-0.44275283,0.03802969,0.00228662},
        {0.00747852,0.01073967,0.02090366,0.00461457,0.03980863,0.00878929,0.00409985,0.00186756,0.18125441,0.00794180,0.01023445,0.00450612,0.01643896,0.26654152,0.00306072,0.01368064,0.00839668,0.10764993,-0.71435091,0.00851526},
        {0.17617706,0.01181629,0.00578676,0.00262530,0.13547871,0.01454379,0.01694332,0.00530363,0.00822937,0.73635171,0.11773937,0.01280613,0.13129028,0.04526924,0.02050210,0.00680190,0.15130413,0.01310401,0.01723920,-1.33539639}
};

AbsTransitionMatrix(inline double)::pythag(double a, double b) {
  double absa = fabs(a), absb = fabs(b);
  return (absa > absb) ?
       absa * sqrt(1+ (absb/absa)*(absb/absa)) :
    absb == 0 ?
       0 :
       absb * sqrt(1+ (absa/absb)*(absa/absb));
}

AbsTransitionMatrix(void)::tred2(double a[], int n, int np, double d[], double e[]){
#define a(i,j) a[(j-1)*np + (i-1)]
#define e(i)   e[i-1]
#define d(i)   d[i-1]
  int i, j, k, l;
  double f, g, h, hh, scale;
  for (i = n; i > 1; i--) {
    l = i-1;
    h = 0;
    scale = 0;
    if ( l > 1 ) {
      for ( k = 1; k <= l; k++ )
	scale += fabs(a(i,k));
      if (scale == 0)
	e(i) = a(i,l);
      else {
	for (k = 1; k <= l; k++) {
	  a(i,k) /= scale;
	  h += a(i,k) * a(i,k);
	}
	f = a(i,l);
	g = -sqrt(h);
	if (f < 0) g = -g;
	e(i) = scale *g;
	h -= f*g;
	a(i,l) = f-g;
	f = 0;
	for (j = 1; j <=l ; j++) {
	  a(j,i) = a(i,j) / h;
	  g = 0;
	  for (k = 1; k <= j; k++)
	    g += a(j,k)*a(i,k);
	  for (k = j+1; k <= l; k++)
	    g += a(k,j)*a(i,k);
	  e(j) = g/h;
	  f += e(j)*a(i,j);
	}
	hh = f/(h+h);
	for (j = 1; j <= l; j++) {
	  f = a(i,j);
	  g = e(j) - hh * f;
	  e(j) = g;
	  for (k = 1; k <= j; k++)
	    a(j,k) -= f*e(k) + g*a(i,k);
	}
      }
    } else
      e(i) = a(i,l);
    d(i) = h;
  }
  d(1) = 0;
  e(1) = 0;
  for (i = 1; i <= n; i++) {
    l = i-1;
    if (d(i) != 0) {
      for (j = 1; j <=l; j++) {
	g = 0;
	for (k = 1; k <= l; k++)
	  g += a(i,k)*a(k,j);
	for (k=1; k <=l; k++)
	  a(k,j) -= g * a(k,i);
      }
    }
    d(i) = a(i,i);
    a(i,i) = 1;
    for (j=1; j<=l; j++)
      a(i,j) = a(j,i) = 0;
  }

  return;
#undef a
#undef e
#undef d
}

AbsTransitionMatrix(void)::tqli(double d[], double e[], int n, int np, double z[]){
#define z(i,j) z[(j-1)*np + (i-1)]
#define e(i)   e[i-1]
#define d(i)   d[i-1]

  int i = 0, iter = 0, k = 0, l = 0, m = 0;
  double b = 0, c = 0, dd = 0, f = 0, g = 0, p = 0, r = 0, s = 0;

  for(i=2; i<=n; i++)
    e(i-1) = e(i);
  e(n) = 0;

  for (l = 1; l <= n; l++)
    {
      iter = 0;
    labelExtra:

      for (m = l; (m < n); m++)
	{
	  dd = fabs(d(m))+fabs(d(m+1));

	  if (fabs(e(m))+dd == dd)
	    break;
	}

      if (m != l)
	{
	  assert(iter < 30);

	  iter++;
	  g = (d(l+1)-d(l))/(2*e(l));
	  r = pythag(g,1.);
	  g = d(m)-d(l)+e(l)/(g+(g<0?-r:r));
	  s = 1;
	  c = 1;
	  p = 0;

	  for (i = m-1; i>=l; i--)
	    {
	      f = s*e(i);
	      b = c*e(i);
	      r = pythag(f,g);

	      e(i+1) = r;
	      if (r == 0)
		{
		  d (i+1) -= p;
		  e (m) = 0;

		  goto labelExtra;
		}
	      s = f/r;
	      c = g/r;
	      g = d(i+1)-p;
	      r = (d(i)-g)*s + 2*c*b;
	      p = s*r;
	      d(i+1) = g + p;
	      g = c*r - b;
	      for (k=1; k <= n; k++)
		{
		  f = z(k,i+1);
		  z(k,i+1) = s * z(k,i) + c*f;
		  z(k,i) = c * z(k,i) - s*f;
		}
	    }
	  d(l) -= p;
	  e(l) = g;
	  e(m) = 0;

	  goto labelExtra;
	}
    }

  return;
#undef z
#undef e
#undef d

}

#undef  double_alignas
#undef AbsTransitionMatrix
