

# generate root dictionary
ROOT_GENERATE_DICTIONARY(G__Dict KKHitInfo.hh NOINSTALL
    LINKDEF LinkDef.h
)
add_library(UnitTests SHARED  G__Dict)
target_include_directories(UnitTests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(UnitTests ${ROOT_LIBRARIES})


# find all unit test sources and build an executable for each
file( GLOB TEST_APP_SOURCES *_unit.cc )

foreach( testsourcefile ${TEST_APP_SOURCES} )
    # get the name of the unit test from the source path
    # (strip _unit.cc, and the path preceding the basename)
    string( REPLACE "_unit.cc" "" testnamenoext ${testsourcefile} )
    get_filename_component(testname ${testnamenoext} NAME)

    # prepend UnitTest_ to the target name to avoid clashes.
    add_executable( UnitTest_${testname} ${testsourcefile} )

    # ensure the unit test filename is just its test name
    set_target_properties( UnitTest_${testname} PROPERTIES OUTPUT_NAME ${testname})

    # link all unit tests to KinKal, MatEnv, and ROOT libraries.
    target_link_libraries( UnitTest_${testname} KinKal MatEnv ${ROOT_LIBRARIES} )

    # register the target as a test 
    add_test (NAME ${testname} COMMAND UnitTest_${testname} )
endforeach( testsourcefile ${TEST_APP_SOURCES} )