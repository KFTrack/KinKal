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
add_library(Examples SHARED ExamplesDict)

target_include_directories(Examples PRIVATE ${PROJECT_SOURCE_DIR}/..)
target_include_directories(Examples PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# link ROOT libraries
target_link_libraries(Examples ${ROOT_LIBRARIES})
