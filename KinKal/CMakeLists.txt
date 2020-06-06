file( GLOB KinKal_LIB_SOURCES *.cc )
file( GLOB KinKal_LIB_HEADERS *.hh )

add_library(KinKal SHARED ${KinKal_LIB_SOURCES} ${KinKal_LIB_HEADERS} )
target_compile_options(KinKal PUBLIC "$<$<CONFIG:DEBUG>:${DBG_FLAGS}>")
target_compile_options(KinKal PUBLIC "$<$<CONFIG:RELEASE>:${PROF_FLAGS}>")
