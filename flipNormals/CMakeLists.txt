
cmake_minimum_required(VERSION 2.8)
 
PROJECT(flipNormals)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(flipNormals MACOSX_BUNDLE flipNormals.cxx )
 
target_link_libraries(flipNormals ${VTK_LIBRARIES})
