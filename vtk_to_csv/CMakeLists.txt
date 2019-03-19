cmake_minimum_required(VERSION 2.8)

set(CMAKE_INSTALL_RPATH "/usr/local/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

PROJECT(vtk_to_csv)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(vtk_to_csv MACOSX_BUNDLE vtk_to_csv.cxx)
 
if(VTK_LIBRARIES)
  target_link_libraries(vtk_to_csv ${VTK_LIBRARIES})
else()
  target_link_libraries(vtk_to_csv vtkHybrid vtkWidgets)
endif()
