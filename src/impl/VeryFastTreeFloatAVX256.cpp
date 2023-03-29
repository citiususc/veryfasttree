#ifdef __AVX__
#include "../Utils.h"
#include "../operations/AVX256Operations.h"
#include "../VeyFastTreeImpl.h"

template class veryfasttree::VeyFastTreeImpl<float, veryfasttree::AVX256Operations>;

#endif

