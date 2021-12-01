# List of unit test sources

set( TEST_SOURCE_FILES
    CentralHelixBField_unit.cc
    CentralHelixDerivs_unit.cc
    CentralHelixFit_unit.cc
    CentralHelixHit_unit.cc
    CentralHelixPKTraj_unit.cc
    CentralHelixTPoca_unit.cc
    CentralHelix_unit.cc
    KinematicLineBField_unit.cc
    KinematicLineDerivs_unit.cc
    KinematicLineFit_unit.cc
    KinematicLineHit_unit.cc
    KinematicLinePKTraj_unit.cc
    KinematicLineTPoca_unit.cc
    KinematicLine_unit.cc
    LoopHelixBField_unit.cc
    LoopHelixDerivs_unit.cc
    LoopHelixFit_unit.cc
    LoopHelixHit_unit.cc
    LoopHelixPKTraj_unit.cc
    LoopHelixTPoca_unit.cc
    LoopHelix_unit.cc
    MatEnv_unit.cc
)

# generate root dictionary
ROOT_GENERATE_DICTIONARY(TestsDict
    # headers to include in ROOT dictionary
    HitInfo.hh
    BFieldInfo.hh
    ParticleTrajectoryInfo.hh
    MaterialInfo.hh

    LINKDEF LinkDef.h
)
# create shared library with ROOT dict
add_library(Tests SHARED TestsDict)

target_include_directories(Tests PRIVATE ${PROJECT_SOURCE_DIR}/..)
target_include_directories(Tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# link ROOT libraries
target_link_libraries(Tests ${ROOT_LIBRARIES})

# Generate unit test targets
foreach( testsourcefile ${TEST_SOURCE_FILES} )
    # get the name of the unit test from the source path
    # (strip _unit.cc, and the path preceding the basename)
    string( REPLACE "_unit.cc" "" testnamenoext ${testsourcefile} )
    get_filename_component(testname ${testnamenoext} NAME)

    # prepend Test_ to the target name to avoid target name clashes,
    # create the unit test executable:
    add_executable( Test_${testname} ${testsourcefile} )

    # add the project root as an include directory
    # link all unit tests to KinKal, MatEnv, and ROOT libraries.
    target_include_directories(Test_${testname} PRIVATE ${PROJECT_SOURCE_DIR}/..)
    target_link_libraries( Test_${testname} General Trajectory Detector Fit MatEnv Tests ${ROOT_LIBRARIES} )

    # ensure the unit test executable filename is just its test name
    set_target_properties( Test_${testname} PROPERTIES OUTPUT_NAME ${testname})

    # register the target as a test
    add_test (NAME ${testname} COMMAND Test_${testname} )
    set_tests_properties(${testname} PROPERTIES TIMEOUT 250)
    set_tests_properties(${testname} PROPERTIES ENVIRONMENT "KINKAL_SOURCE_DIR=${CMAKE_SOURCE_DIR}")

    install( TARGETS Test_${testname}
             RUNTIME DESTINATION bin/ )

endforeach( testsourcefile ${TEST_SOURCE_FILES} )
