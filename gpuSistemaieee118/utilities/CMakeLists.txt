include_directories(${CUDA_SISTEMA_IEEE_FC_118_SOURCE_DIR}/utilities)
link_directories(${CUDA_SISTEMA_IEEE_FC_118_BINARY_DIR}/utilities)

cuda_add_library(utilities utilities.cu)
cuda_add_library(cusolverUtilities cusolverUtilities.cu)
