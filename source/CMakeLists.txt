CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(analytical_composite_regions)

INCLUDE_DIRECTORIES($ENV{GITHOME}/libraries)

ADD_LIBRARY(mylib $ENV{GITHOME}/libraries/DataTools.cpp $ENV{GITHOME}/libraries/Names.cpp)

ADD_EXECUTABLE(mycode composite_regions_6_6.cc)

SET(LAPACK_DIR "" CACHE PATH "path to liblapack")

TARGET_LINK_LIBRARIES(mycode mylib lapack blas)

TARGET_LINK_DIRECTORIES(mycode PUBLIC ${LAPACK_DIR})
