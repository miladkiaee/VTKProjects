
cmake_minimum_required(VERSION 2.8)
 
PROJECT(PolyDataExtractNormals)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(PolyDataExtractNormals MACOSX_BUNDLE PolyDataExtractNormals.cxx )
 
target_link_libraries(PolyDataExtractNormals ${VTK_LIBRARIES})
