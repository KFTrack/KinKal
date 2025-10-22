
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
)

# set top-level directory as include root
target_include_directories(MatEnv PUBLIC 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:include>)

# set shared library version equal to project version
set_target_properties(MatEnv PROPERTIES VERSION ${PROJECT_VERSION} PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
