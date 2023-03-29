#ifdef __AVX512F__
#include "../Utils.h"
#include "../operations/AVX512Operations.h"
#include "../VeyFastTreeImpl.h"

template class veryfasttree::VeyFastTreeImpl<double, veryfasttree::AVX512Operations>;

#endif