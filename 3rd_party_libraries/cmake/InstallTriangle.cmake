# triangle package

include(ExternalProject)

set(TRIANGLE_SOURCE_DIR ${MAIN_SOURCE_DIR}/triangle)
set(TRIANGLE_BINARY_DIR ${MAIN_BINARY_DIR}/triangle)
set(TRIANGLE_INSTALL_PREFIX ${MAIN_INSTALL_PREFIX}/triangle)

message(STATUS "Install triangle release 1.0.4")

ExternalProject_Add(triangle
    GIT_REPOSITORY https://github.com/fvicini/triangle.git
    GIT_TAG ecda41addc40b79429e491924dafa80a00887b99 # release 1.0.4
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    SOURCE_DIR ${TRIANGLE_SOURCE_DIR}
    BINARY_DIR ${TRIANGLE_BINARY_DIR}
    CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${TRIANGLE_INSTALL_PREFIX}
    )
