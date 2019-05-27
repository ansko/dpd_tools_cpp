cmake_minimum_required (VERSION 2.6) 
 

set (PROJECT dpd_tools) 


project (${PROJECT}) 


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
#
#   -DMACRO define MACRO:
#     DEBUG               output debug information
#     DPARTIAL_DATAFILES  write out some incomplete datafiles (for debug mainly)

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}\
            -Werror \
            -std=c++17 \
            -g \
            -lm \
            -O3 \
            -DDEBUG \
    ")
endif (CMAKE_COMPILER_IS_GNUCXX)


set (HEADERS
    src/structure.hpp
    src/options.hpp
    src/write_data.hpp
)


add_executable (isolated ${HEADERS} ${SOURCES} isolated.cpp)
add_executable (isolated_anisotropic ${HEADERS} ${SOURCES} isolated_anisotropic.cpp)
add_executable (periodic ${HEADERS} ${SOURCES} periodic.cpp)
