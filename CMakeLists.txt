cmake_minimum_required(VERSION 3.9)
project(GeDiMProject LANGUAGES CXX C VERSION 1.0)

# Include library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/GeDiM)

# Include unit tests
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/GeDiM/test)
