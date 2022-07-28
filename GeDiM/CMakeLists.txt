cmake_minimum_required(VERSION 3.9)

project(GeDiM LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
###########################################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create library variables
###############################################################################################
set(GEDIM_LIBRARY_SOURCES ${PROJECT_NAME}_sources)
set(GEDIM_LIBRARY_HEADERS ${PROJECT_NAME}_headers)
set(GEDIM_LIBRARY_PUBLIC_HEADERS ${PROJECT_NAME}_public_headers)
set(GEDIM_LIBRARY_LINKED_LIBRARIES ${PROJECT_NAME}_linked_libraries)
set(GEDIM_LIBRARY_INCLUDE ${PROJECT_NAME}_include)
set(GEDIM_LIBRARY_ADDITIONAL_INCLUDE ${PROJECT_NAME}_additional_include)
set(GEDIM_INSTALL_DIR ${PROJECT_NAME})
set(GEDIM_INSTALL_BINARY_DIR bin/)
set(GEDIM_INSTALL_INCLUDE_DIR include/)
set(GEDIM_INSTALL_LIB_DIR lib/)
set(GEDIM_INSTALL_ARCHIVE_DIR lib/)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
Message(STATUS "Gedim Build configuration: " ${CMAKE_BUILD_TYPE})

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${GEDIM_INSTALL_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
Message(STATUS "Gedim Library will be installed in: " ${CMAKE_INSTALL_PREFIX})

set(USE_MPI OFF CACHE BOOL "Enable MPI support")
set(VERBOSE 3 CACHE STRING "Select verbosity level (0, 1, 2 3 or 4)")
set_property(CACHE VERBOSE PROPERTY STRINGS 0 1 2 3 4)
if(NOT ${VERBOSE} MATCHES 0 AND NOT ${VERBOSE} MATCHES 1 AND NOT ${VERBOSE} MATCHES 2 AND NOT ${VERBOSE} MATCHES 3 AND NOT ${VERBOSE} MATCHES 4)
        message(SEND_ERROR "VERBOSE variable must be either 0 (None), 1 (Error), 2 (Warning), 3 (Info) or 4 (Debug)")
endif()

# Add dependencies
###########################################################################################
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND ${GEDIM_LIBRARY_LINKED_LIBRARIES} PUBLIC Eigen3::Eigen)

# Insert Sources
###############################################################################################
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
list(APPEND ${GEDIM_LIBRARY_SOURCES} ${gedim_src_sources})
list(APPEND ${GEDIM_LIBRARY_HEADERS} ${gedim_src_headers})
list(APPEND ${GEDIM_LIBRARY_INCLUDE} ${gedim_src_includes})

# Create public headers
###############################################################################################
set(GEDIM_LIBRARY_ADDITIONAL_INCLUDE_PATH ${CMAKE_CURRENT_BINARY_DIR}/${GEDIM_LIBRARY_ADDITIONAL_INCLUDE})

foreach(header ${${GEDIM_LIBRARY_HEADERS}})
        get_filename_component(headerDirectory ${header} DIRECTORY)
        get_filename_component(headerFileName ${header} NAME_WE)
        get_filename_component(headerExtension ${header} EXT)

  if(${headerExtension} MATCHES ".hpp.in")
    set(headerExtension ".hpp")
    set(headerDirectory ${GEDIM_LIBRARY_ADDITIONAL_INCLUDE_PATH})

    configure_file(${header} ${headerDirectory}/${headerFileName}${headerExtension} @ONLY)
  elseif(NOT ${headerExtension} MATCHES ".hpp" AND NOT ${headerExtension} MATCHES ".h")
    message(SEND_ERROR "Header extension not recognized for file ${headerFileName}${headerExtension}")
  endif()

  list(APPEND ${GEDIM_LIBRARY_PUBLIC_HEADERS} ${headerDirectory}/${headerFileName}${headerExtension})
endforeach()

list(APPEND ${GEDIM_LIBRARY_INCLUDE} ${GEDIM_LIBRARY_ADDITIONAL_INCLUDE_PATH})

# Create library
###############################################################################################
add_library(${PROJECT_NAME} STATIC ${${GEDIM_LIBRARY_SOURCES}} ${${GEDIM_LIBRARY_HEADERS}})
add_library(GeDiM::GeDiM ALIAS ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} ${${GEDIM_LIBRARY_LINKED_LIBRARIES}})

target_include_directories(${PROJECT_NAME} PRIVATE ${${GEDIM_LIBRARY_INCLUDE}})
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:${GEDIM_INSTALL_INCLUDE_DIR}>)

target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

# Create install target
###############################################################################################
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${${GEDIM_LIBRARY_PUBLIC_HEADERS}}")
set_target_properties(${PROJECT_NAME} PROPERTIES EXPORT_NAME GeDiM)
install(
	TARGETS ${PROJECT_NAME} 
	EXPORT GeDiMTargets
    RUNTIME DESTINATION ${GEDIM_INSTALL_BINARY_DIR} # runtime executable
    LIBRARY DESTINATION ${GEDIM_INSTALL_LIB_DIR} # dynamic libraries
    ARCHIVE DESTINATION ${GEDIM_INSTALL_ARCHIVE_DIR} # static libraries
    PUBLIC_HEADER DESTINATION ${GEDIM_INSTALL_INCLUDE_DIR} # headers
)


include(CMakePackageConfigHelpers)

configure_package_config_file(
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/GeDiMConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/GeDiMConfig.cmake
	INSTALL_DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${GEDIM_INSTALL_DIR}/cmake
	PATH_VARS GEDIM_INSTALL_INCLUDE_DIR
  )

write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/GeDiMConfigVersion.cmake
	VERSION ${GeDiM_VERSION}
	COMPATIBILITY AnyNewerVersion
)

export(
	TARGETS ${PROJECT_NAME} 
	NAMESPACE GeDiM::
	FILE ${CMAKE_CURRENT_BINARY_DIR}/GeDiMTargets.cmake
)

install(
	EXPORT GeDiMTargets 
	NAMESPACE GeDiM:: 
	DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${GEDIM_INSTALL_DIR}/cmake
)

install(
	FILES
		${CMAKE_CURRENT_BINARY_DIR}/GeDiMConfig.cmake
		${CMAKE_CURRENT_BINARY_DIR}/GeDiMConfigVersion.cmake
	DESTINATION
		${CMAKE_CURRENT_BINARY_DIR}/${GEDIM_INSTALL_DIR}/cmake
)