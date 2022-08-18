# VTK package

include(ExternalProject)

set(VTK_SOURCE_DIR ${MAIN_SOURCE_DIR}/vtk)
set(VTK_BINARY_DIR ${MAIN_BINARY_DIR}/vtk)
set(VTK_INSTALL_PREFIX ${MAIN_INSTALL_PREFIX}/vtk)

message(STATUS "Install VTK release 9.1.0")

ExternalProject_Add(VTK
    GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git
    GIT_TAG 285daeedd58eb890cb90d6e907d822eea3d2d092 # release 9.1.0
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    SOURCE_DIR ${VTK_SOURCE_DIR}
    BINARY_DIR ${VTK_BINARY_DIR}
    CMAKE_ARGS -DVTK_GROUP_ENABLE_Imaging=NO -DVTK_GROUP_ENABLE_MPI=NO -DVTK_GROUP_ENABLE_Qt=NO -DVTK_GROUP_ENABLE_StandAlone=NO -DVTK_GROUP_ENABLE_Views=NO -DVTK_GROUP_ENABLE_Web=NO -DVTK_GROUP_ENABLE_StandAlone=NO -DBUILD_SHARED_LIBS=ON -DVTK_MODULE_ENABLE_VTK_CommonCore=YES -DVTK_MODULE_ENABLE_VTK_CommonExecutionModel=YES -DVTK_MODULE_ENABLE_VTK_IOXML=YES -DVTK_MODULE_ENABLE_VTK_IOXMLParser=YES -DVTK_MODULE_ENABLE_VTK_CommonMisc=YES -DVTK_MODULE_ENABLE_VTK_CommonSystem=YES -DVTK_MODULE_ENABLE_VTK_IOCore=YES -DVTK_MODULE_ENABLE_VTK_CommonMath=YES -DVTK_MODULE_ENABLE_VTK_CommonTransforms=YES -DVTK_MODULE_ENABLE_VTK_FiltersCore=YES -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${VTK_INSTALL_PREFIX}
    )
