#if (defined __SSE2__) || (defined __AVX__)
#include "../Utils.h"
#include "../operations/SSE128Operations.h"
#include "../VeyFastTreeImpl.h"

template class veryfasttree::VeyFastTreeImpl<float, veryfasttree::SSE128Operations>;

#endif