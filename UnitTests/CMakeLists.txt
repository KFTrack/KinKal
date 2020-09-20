# List of unit test sources
# ls -1 *_unit.cc

set( TEST_SOURCE_FILES
    CentralHelixBFieldTest_unit.cc
    CentralHelixDerivs_unit.cc
    CentralHelixFitTest_unit.cc
    CentralHelixHitTest_unit.cc
    CentralHelixPKTrajTest_unit.cc
    CentralHelixTPocaTest_unit.cc
    CentralHelix_unit.cc
    KinematicLineBFieldTest_unit.cc
    KinematicLineDerivs_unit.cc
    KinematicLineFitTest_unit.cc
    KinematicLineHitTest_unit.cc
    KinematicLinePKTrajTest_unit.cc
    KinematicLineTPocaTest_unit.cc
    KinematicLine_unit.cc
    LoopHelixBFieldTest_unit.cc
    LoopHelixDerivs_unit.cc
    LoopHelixFitTest_unit.cc
    LoopHelixHitTest_unit.cc
    LoopHelixPKTrajTest_unit.cc
    LoopHelixTPocaTest_unit.cc
    LoopHelix_unit.cc
    MatEnv_unit.cc
)



# generate root dictionary

ROOT_GENERATE_DICTIONARY(G__Dict 
    # headers to include in ROOT dictionary
    HitInfo.hh 
    BFieldInfo.hh 
    ParticleTrajectoryInfo.hh 
    MaterialInfo.hh 
    NOINSTALL
    LINKDEF LinkDef.h
)
# create shared library with ROOT dict
add_library(UnitTests SHARED G__Dict)

target_include_directories(UnitTests PRIVATE ${PROJECT_SOURCE_DIR})
target_include_directories(UnitTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# link ROOT libraries
target_link_libraries(UnitTests ${ROOT_LIBRARIES})

# Generate unit test targets
foreach( testsourcefile ${TEST_SOURCE_FILES} )
    # get the name of the unit test from the source path
    # (strip _unit.cc, and the path preceding the basename)
    string( REPLACE "_unit.cc" "" testnamenoext ${testsourcefile} )
    get_filename_component(testname ${testnamenoext} NAME)

    # prepend UnitTest_ to the target name to avoid target name clashes,
    # create the unit test executable:
    add_executable( UnitTest_${testname} ${testsourcefile} )

    # add the project root as an include directory
    # link all unit tests to KinKal, MatEnv, and ROOT libraries.
    target_include_directories(UnitTest_${testname} PRIVATE ${PROJECT_SOURCE_DIR})
    target_link_libraries( UnitTest_${testname} General Trajectory Detector Fit MatEnv UnitTests ${ROOT_LIBRARIES} )

    # ensure the unit test executable filename is just its test name
    set_target_properties( UnitTest_${testname} PROPERTIES OUTPUT_NAME ${testname})

    # register the target as a test 
    add_test (NAME ${testname} COMMAND UnitTest_${testname} )
    set_tests_properties(${testname} PROPERTIES TIMEOUT 15) 
    set_tests_properties(${testname} PROPERTIES ENVIRONMENT "PACKAGE_SOURCE=${CMAKE_SOURCE_DIR}")

    install( TARGETS UnitTest_${testname}
             RUNTIME DESTINATION bin/ )
 
endforeach( testsourcefile ${TEST_SOURCE_FILES} )
