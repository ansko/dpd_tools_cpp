cmake_minimum_required (VERSION 2.6) 
 

set (PROJECT dpd_tools) 


project (${PROJECT}) 


if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -std=c++17 -g -lm  -O3 -DDEBUG")
endif (CMAKE_COMPILER_IS_GNUCXX)


set (HEADERS
    src/structure.hpp
    src/options.hpp
    src/write_data.hpp
)
set (SOURCES
    #src/tactoid.cpp
)


add_executable (tactoid ${HEADERS} ${SOURCES} tactoid.cpp)
