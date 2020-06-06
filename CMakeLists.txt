cmake_minimum_required(VERSION 2.8.12)

project (KinKal)

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