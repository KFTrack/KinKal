# Explicitly list source files in this shared library
# When a new source file is added, it must be added to this list.

# generate root dictionary
ROOT_GENERATE_DICTIONARY(ExamplesDict
    # headers to include in ROOT dictionary
    HitInfo.hh
    BFieldInfo.hh
    ParticleTrajectoryInfo.hh
    MaterialInfo.hh
    LINKDEF LinkDef.h
    MODULE Examples
)
# create shared library with ROOT dict
add_library(Examples SHARED
  CaloDistanceToTime.cc
  ExamplesDict
)
target_include_directories(Examples PRIVATE ${PROJECT_SOURCE_DIR}/..)
target_include_directories(Examples PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# link ROOT libraries
target_link_libraries(Examples ${ROOT_LIBRARIES})
# set shared library version equal to project version
set_target_properties(Examples PROPERTIES VERSION ${PROJECT_VERSION} PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
