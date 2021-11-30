
add_library(MatEnv SHARED
    DetMaterial.cc
    ElmPropObj.cc
    MatDBInfo.cc
    MatElementList.cc
    MatElementObj.cc
    MatElmDictionary.cc
    MatIsoDictionary.cc
    MatIsotopeList.cc
    MatIsotopeObj.cc
    MatMaterialList.cc
    MatMaterialObj.cc
    MatMtrDictionary.cc
    MtrPropObj.cc
    RecoMatFactory.cc
    SimpleFileFinder.cc
    MoyalDist.cc
)

# set top-level directory as include root
target_include_directories(MatEnv PRIVATE ${PROJECT_SOURCE_DIR}/..)

# set shared library version equal to project version
set_target_properties(MatEnv PROPERTIES VERSION ${PROJECT_VERSION})
