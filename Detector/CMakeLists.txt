
# Explicitly list source files in this shared library
# When a new source file is added, it must be added to this list.

# Globbing is not recommended, see: https://cmake.org/cmake/help/v3.12/command/file.html#filesystem

# you can regenerate this list easily by running in this directory: ls -1 *.cc
add_library(Detector SHARED
    StrawMaterial.cc
    Residual.cc
)

#message( "source dir detector " ${CMAKE_SOURCE_DIR}/..)
# set top-level directory as include root
target_include_directories(Detector PRIVATE ${PROJECT_SOURCE_DIR}/..)

# link this library with ROOT libraries and other KinKal libraries
target_link_libraries(Detector Trajectory General ${ROOT_LIBRARIES})

# set shared library version equal to project version
set_target_properties(Detector PROPERTIES VERSION ${PROJECT_VERSION})
