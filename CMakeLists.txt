cmake_minimum_required(VERSION 3.9.0)

# always use gcc
set(CMAKE_C_COMPILER "gcc")

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Og")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message( "Build Type: ${CMAKE_BUILD_TYPE}" )

project (KinKal CXX)

find_package(ROOT REQUIRED COMPONENTS Core RIO Net Hist Spectrum GenVector MLP Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread Gui)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(MatEnv)
add_subdirectory(KinKal)
add_subdirectory(UnitTests)