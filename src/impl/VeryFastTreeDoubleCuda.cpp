#if USE_CUDA
#include "../Utils.h"
#include "../operations/CudaOperations.h"
#include "../VeyFastTreeImpl.h"

template class veryfasttree::VeyFastTreeImpl<double, veryfasttree::CudaOperations>;

#endif