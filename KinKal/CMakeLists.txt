
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

# select which headers to install ls -1 *.hh
install(FILES
    BFieldEffect.hh
    BFieldMap.hh
    BFieldUtils.hh
    CentralHelix.hh
    ClosestApproachData.hh
    ClosestApproach.hh
    Config.hh
    Constraint.hh
    DetectorHit.hh
    DetectorXing.hh
    Effect.hh
    FitData.hh
    FitState.hh
    FitStatus.hh
    KinematicLine.hh
    Line.hh
    LoopHelix.hh
    LRAmbig.hh
    Material.hh
    MaterialXing.hh
    Measurement.hh
    MomBasis.hh
    Parameters.hh
    ParticleState.hh
    ParticleTrajectory.hh
    PieceClosestApproach.hh
    PieceTrajectory.hh
    POCAUtil.hh
    PrintDetail.hh
    Residual.hh
    ScintHit.hh
    StrawHit.hh
    StrawMat.hh
    StrawXing.hh
    TimeDir.hh
    TimeRange.hh
    TrackEnd.hh
    Track.hh
    Vectors.hh
    Weights.hh
    WireCell.hh
    WireHit.hh

    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/KinKal
)
