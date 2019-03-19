
cmake_minimum_required(VERSION 2.8)
 
PROJECT(WindowedSincPolyDataFilter)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(WindowedSincPolyDataFilter MACOSX_BUNDLE WindowedSincPolyDataFilter.cxx )
 
target_link_libraries(WindowedSincPolyDataFilter ${VTK_LIBRARIES})
