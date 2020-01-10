CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(analytical_composite_regions)

INCLUDE_DIRECTORIES(/home/$ENV{USER}/git/libraries)

ADD_LIBRARY(mylib /home/$ENV{USER}/git/libraries/DataTools.cpp /home/$ENV{USER}/git/libraries/Names.cpp)

ADD_EXECUTABLE(mycode composite_regions_6_6.cc)

TARGET_LINK_LIBRARIES(mycode mylib lapack blas)

