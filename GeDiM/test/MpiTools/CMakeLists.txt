list(APPEND gedim_test_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_metisUtilities.hpp)
list(APPEND gedim_test_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_graphUtilities.hpp)

list(APPEND gedim_test_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(gedim_test_src_sources ${gedim_test_src_sources} PARENT_SCOPE)
set(gedim_test_src_headers ${gedim_test_src_headers} PARENT_SCOPE)
set(gedim_test_src_includes ${gedim_test_src_includes} PARENT_SCOPE)
