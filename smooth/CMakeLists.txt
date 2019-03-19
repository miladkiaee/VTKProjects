
cmake_minimum_required(VERSION 2.8)
 
PROJECT(smooth)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(smooth MACOSX_BUNDLE smooth.cxx )
 
target_link_libraries(smooth ${VTK_LIBRARIES})
