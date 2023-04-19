#if USE_CUDA
#include "../Utils.h"
#include "../operations/CudaOperations.h"
#include "../VeyFastTreeImpl.h"

template class veryfasttree::VeyFastTreeImpl<float, veryfasttree::CudaOperations>;

#endif