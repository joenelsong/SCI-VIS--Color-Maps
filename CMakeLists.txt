cmake_minimum_required(VERSION 2.8)

PROJECT(proj3)
SET(VTK_DIR /home/operatingsystems/VTK-6.3.0)

find_package(VTK REQUIRED)
#find_package(Qt5Widgets)
include(${VTK_USE_FILE})

add_executable(proj3 proj3)

if(VTK_LIBRARIES)
  target_link_libraries(proj3 ${VTK_LIBRARIES})
else()
  target_link_libraries(proj3 vtkHybrid)
endif()


