cmake_minimum_required (VERSION 2.6) 
 

set (PROJECT dpd_tools) 


project (${PROJECT}) 

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)


# CMAKE_CXX_FLAGS:
#
#   -std=c++17
#
#   -g0    no debug information
#   -g1    minimal debug information
#   -g     default debug information
#   -g3    maximal debug information
#
#   -lm    link math library

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}\
            -Werror \
            -std=c++17 \
            -g \
            -lm \
            -O3 \
    ")
endif (CMAKE_COMPILER_IS_GNUCXX)

#https://stackoverflow.com/a/9328525/7927226
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()



set (HEADERS
    src/options_parser.hpp
    src/structures.hpp
    src/write_data.hpp
    src/parameters_classes.hpp
)

set (SOURCES
    src/structure_add_mmt_circular.cpp
    src/structure_add_mmt_circular_fcc.cpp

    src/structure_add_mmt_periodic.cpp
    src/structure_add_modifier_gallery.cpp
    src/structure_add_poly_bbox.cpp
)


add_executable (isolated_parallel ${HEADERS} ${SOURCES} isolated_parallel.cpp)
target_link_libraries(isolated_parallel Threads::Threads)

add_executable (isolated_parallel_fcc ${HEADERS} ${SOURCES} isolated_parallel_fcc.cpp)
target_link_libraries(isolated_parallel_fcc Threads::Threads)

add_executable (periodic ${HEADERS} ${SOURCES} periodic.cpp)
