
# Explicitly list source files in this shared library
# When a new source file is added, it must be added to this list.

# Globbing is not recommended, see: https://cmake.org/cmake/help/v3.12/command/file.html#filesystem

# you can regenerate this list easily by running in this directory: ls -1 *.cc
add_library(KinKal SHARED 
    BFieldMap.cc
    CentralHelix.cc
    ClosestApproachData.cc
    Config.cc
    FitStatus.cc
    KinematicLine.cc
    Line.cc
    LoopHelix.cc
    LRAmbig.cc
    Parameters.cc
    ParticleState.cc
    POCAUtil.cc
    StrawMat.cc
    TimeDir.cc
    TimeRange.cc
    Weights.cc
)

# set top-level directory as include root
target_include_directories(KinKal PRIVATE ${PROJECT_SOURCE_DIR})

# link this library with ROOT libraries
target_link_libraries(KinKal ${ROOT_LIBRARIES})

# set shared library version equal to project version
set_target_properties(KinKal PROPERTIES VERSION ${PROJECT_VERSION})
