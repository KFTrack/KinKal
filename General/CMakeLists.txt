
# Explicitly list source files in this shared library
# When a new source file is added, it must be added to this list.

add_library(General SHARED
    Chisq.cc
    Parameters.cc
    ParticleState.cc
    TimeDir.cc
    TimeRange.cc
    Weights.cc
    BFieldMap.cc
    AxialBFieldMap.cc
    )

# create shared library with ROOT dict
ROOT_GENERATE_DICTIONARY(GeneralDict
    ParticleState.hh
    Vectors.hh
    LINKDEF LinkDef.h
    MODULE General
)

# set top-level directory as include root
target_include_directories(General PRIVATE ${PROJECT_SOURCE_DIR}/..)

# link this library with ROOT libraries
target_link_libraries(General ${ROOT_LIBRARIES})

# set shared library version equal to project version
set_target_properties(General PROPERTIES VERSION ${PROJECT_VERSION})
