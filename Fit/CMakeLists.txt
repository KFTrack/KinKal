
# Explicitly list source files in this shared library
# When a new source file is added, it must be added to this list.

# Globbing is not recommended, see: https://cmake.org/cmake/help/v3.12/command/file.html#filesystem

# you can regenerate this list easily by running in this directory: ls -1 *.cc
add_library(Fit SHARED
    Config.cc
    MetaIterConfig.cc
    Status.cc
)

# set top-level directory as include root
target_include_directories(Fit PUBLIC 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/..>
    $<INSTALL_INTERFACE:include>)

# link this library with ROOT libraries
target_link_libraries(Fit Detector Trajectory General ${ROOT_LIBRARIES})

# set shared library version equal to project version
set_target_properties(Fit PROPERTIES VERSION ${PROJECT_VERSION} PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
