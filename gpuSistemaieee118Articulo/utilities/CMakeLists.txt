include_directories(${CUDA_SISTEMA_IEEE_FC_118_ARTICULO_SOURCE_DIR}/utilities)
link_directories(${CUDA_SISTEMA_IEEE_FC_118_ARTICULO_BINARY_DIR}/utilities)

cuda_add_library(utilities utilities.c)
cuda_add_library(cusolverUtilities cusolverUtilities.cu)
