
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
)

# set top-level directory as include root
target_include_directories(MatEnv PRIVATE ${PROJECT_SOURCE_DIR})

# set shared library version equal to project version
set_target_properties(MatEnv PROPERTIES VERSION ${PROJECT_VERSION})

# select which headers to install ls -1 *.hh
install(FILES
    BbrCollectionUtils.hh
    DetMaterial.hh
    ElmPropObj.hh
    ErrLog.hh
    FileFinderInterface.hh
    MatDBInfo.hh
    MatElementList.hh
    MatElementObj.hh
    MatElmDictionary.hh
    MaterialInfo.hh
    MatIsoDictionary.hh
    MatIsotopeList.hh
    MatIsotopeObj.hh
    MatMaterialList.hh
    MatMaterialObj.hh
    MatMtrDictionary.hh
    MtrPropObj.hh
    PtrLess.hh
    RecoMatFactory.hh

    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MatEnv
)
