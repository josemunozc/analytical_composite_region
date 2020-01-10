CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(analytical_composite_regions)

INCLUDE_DIRECTORIES($ENV{HOME}/git/libraries)

ADD_LIBRARY(mylib $ENV{HOME}/git/libraries/DataTools.cpp $ENV{HOME}/git/libraries/Names.cpp)

ADD_EXECUTABLE(mycode composite_regions_6_6.cc)

TARGET_LINK_LIBRARIES(mycode mylib lapack blas)

