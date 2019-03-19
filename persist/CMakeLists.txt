cmake_minimum_required(VERSION 2.8)

set(CMAKE_INSTALL_RPATH "/usr/local/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

PROJECT(persist)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(persist MACOSX_BUNDLE persist.cxx)
 
if(VTK_LIBRARIES)
  target_link_libraries(persist ${VTK_LIBRARIES})
else()
  target_link_libraries(persist vtkHybrid vtkWidgets)
endif()
