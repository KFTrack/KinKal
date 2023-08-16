# Explicitly list source files in this shared library
# When a new source file is added, it must be added to this list.

# Globbing is not recommended, see: https://cmake.org/cmake/help/v3.12/command/file.html#filesystem

# you can regenerate this list easily by running in this directory: ls -1 *.cc
add_library(Geometry SHARED
  Annulus.cc
  Cylinder.cc
  Disk.cc
  IntersectFlag.cc
  Plane.cc
  Ray.cc
  Rectangle.cc
)
# make the library names unique
set(CMAKE_SHARED_LIBRARY_PREFIX "libKinKal_")

# set top-level directory as include root
target_include_directories(Geometry PRIVATE ${PROJECT_SOURCE_DIR}/..)

# link this library with ROOT libraries
target_link_libraries(Geometry General ${ROOT_LIBRARIES})

# set shared library version equal to project version
set_target_properties(Geometry PROPERTIES VERSION ${PROJECT_VERSION} PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
