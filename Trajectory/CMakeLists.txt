
# Explicitly list source files in this shared library
# When a new source file is added, it must be added to this list.

# Globbing is not recommended, see: https://cmake.org/cmake/help/v3.12/command/file.html#filesystem

# you can regenerate this list easily by running in this directory: ls -1 *.cc
add_library(Trajectory SHARED 
    CentralHelix.cc
    ClosestApproachData.cc
    KinematicLine.cc
    Line.cc
    LoopHelix.cc
    POCAUtil.cc
)

# set top-level directory as include root
target_include_directories(Trajectory PRIVATE ${PROJECT_SOURCE_DIR}/..)

# link this library with ROOT libraries
target_link_libraries(Trajectory General ${ROOT_LIBRARIES})

# set shared library version equal to project version
set_target_properties(Trajectory PROPERTIES VERSION ${PROJECT_VERSION})
