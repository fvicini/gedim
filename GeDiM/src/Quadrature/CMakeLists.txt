list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/QuadratureData.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_Gauss1D.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_GaussLobatto1D.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_Gauss2D_Triangle.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_Gauss2D_Square.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_Gauss3D_Tetrahedron.hpp)
list(APPEND gedim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_Gauss3D_Hexahedron.hpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_Gauss3D_Tetrahedron_PositiveWeights.hpp)

list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_Gauss2D_Triangle.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_Gauss2D_Square.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_Gauss3D_Tetrahedron.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_Gauss3D_Hexahedron.cpp)
list(APPEND gedim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Quadrature_Gauss3D_Tetrahedron_PositiveWeights.cpp)

list(APPEND gedim_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(gedim_src_sources ${gedim_src_sources} PARENT_SCOPE)
set(gedim_src_headers ${gedim_src_headers} PARENT_SCOPE)
set(gedim_src_includes ${gedim_src_includes} PARENT_SCOPE)
