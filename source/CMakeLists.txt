CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(analytical_composite_regions)

INCLUDE_DIRECTORIES($ENV{GITHOME}/libraries)

ADD_LIBRARY(mylib $ENV{GITHOME}/libraries/DataTools.cpp $ENV{GITHOME}/libraries/Names.cpp)

SET(LAPACK_DIR "" CACHE PATH "path to liblapack")

LINK_DIRECTORIES(${LAPACK_DIR})
ADD_EXECUTABLE(mycode composite_regions_6_6.cc)
TARGET_LINK_LIBRARIES(mycode mylib mylapack blas)

