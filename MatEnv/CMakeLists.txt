
file( GLOB MATENV_LIB_SOURCES *.cc )
file( GLOB MATENV_LIB_HEADERS *.hh )

add_library(MatEnv SHARED ${MATENV_LIB_SOURCES} ${MATENV_LIB_HEADERS} )
target_compile_options(MatEnv PUBLIC "$<$<CONFIG:DEBUG>:${DBG_FLAGS}>")
target_compile_options(MatEnv PUBLIC "$<$<CONFIG:RELEASE>:${PROF_FLAGS}>")
