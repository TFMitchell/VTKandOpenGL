cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project1C)
SET(VTK_DIR c:/VTKbin)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(project1C project1C)

if(VTK_LIBRARIES)
  target_link_libraries(project1C ${VTK_LIBRARIES})
else()
  target_link_libraries(project1C vtkHybrid)
endif()

