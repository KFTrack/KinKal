file( GLOB KinKal_LIB_SOURCES *.cc )
file( GLOB KinKal_LIB_HEADERS *.hh )

add_library(KinKal SHARED ${KinKal_LIB_SOURCES} ${KinKal_LIB_HEADERS} )
