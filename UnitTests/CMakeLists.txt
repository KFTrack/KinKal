
# generate root dictionary

ROOT_GENERATE_DICTIONARY(G__Dict KKHitInfo.hh NOINSTALL
    LINKDEF LinkDef.h
)
add_library(UnitTests SHARED  G__Dict)
target_include_directories(UnitTests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(UnitTests ${ROOT_LIBRARIES})

file( GLOB TEST_APP_SOURCES *_unit.cc )

foreach( testsourcefile ${TEST_APP_SOURCES} )
    string( REPLACE ".cc" "" testnamenoext ${testsourcefile} )
    get_filename_component(testname ${testnamenoext} NAME)

    add_executable( ${testname} ${testsourcefile} )

    target_compile_options(${testname} PUBLIC "$<$<CONFIG:DEBUG>:${DBG_FLAGS}>")
    target_compile_options(${testname} PUBLIC "$<$<CONFIG:RELEASE>:${PROF_FLAGS}>")

    target_link_libraries( ${testname} KinKal MatEnv ${ROOT_LIBRARIES} )
endforeach( testsourcefile ${TEST_APP_SOURCES} )