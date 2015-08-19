include_directories(${SISTEMA_IEEE_FC_118_ARTICULO_SOURCE_DIR}/utilities)
link_directories(${SISTEMA_IEEE_FC_118_ARTICULO_BINARY_DIR}/utilities)

add_library(utilities utilities.c)
