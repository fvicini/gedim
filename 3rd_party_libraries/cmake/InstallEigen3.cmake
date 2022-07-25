# Eigen3 package

include(ExternalProject)

set(EIGEN3_SOURCE_DIR ${MAIN_SOURCE_DIR}/eigen3)
set(EIGEN3_BINARY_DIR ${MAIN_BINARY_DIR}/eigen3)
set(EIGEN3_INSTALL_PREFIX ${MAIN_INSTALL_PREFIX}/eigen3)

message(STATUS "Install Eigen3 release 3.4.0")

ExternalProject_Add(Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3147391d946bb4b6c68edd901f2add6ac1f31f8c # release 3.4.0
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    SOURCE_DIR ${EIGEN3_SOURCE_DIR}
    BINARY_DIR ${EIGEN3_BINARY_DIR}
    CMAKE_ARGS -DBUILD_TESTING=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_PKGCONFIG=OFF -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${EIGEN3_INSTALL_PREFIX}
    )
