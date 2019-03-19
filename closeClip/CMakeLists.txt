
cmake_minimum_required(VERSION 2.8)
 
PROJECT(closeClip)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(closeClip MACOSX_BUNDLE closeClip.cxx )
 
target_link_libraries(closeClip ${VTK_LIBRARIES})
