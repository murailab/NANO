cmake_minimum_required(VERSION 3.13)
project(nano VERSION 0.2
        DESCRIPTION "Analysis of Astrocytic nanoarchitecture"
        LANGUAGES CXX)	
message(STATUS "Found CMAKE-${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}")
message(STATUS "----------------Running CMake---------------------------------")

#set(ITK_DIR <PATH_OF_ITK_INSTALL>/itk/lib/cmake/ITK-5.X)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-g -rdynamic -O3 -Wall")

message(STATUS "CMAKE_CXX_FLAGS : ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_STANDARD : ${CMAKE_CXX_STANDARD}")

add_library(daedalum)

add_subdirectory(common)
add_subdirectory(marching)

add_executable(nano main.cpp)

if((${CMAKE_COMPILER_IS_GNUCXX}) AND (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS_EQUAL 8))
        message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}: adding stdc++fs library to link")
        target_link_libraries(nano daedalum ${ITK_LIBRARIES} stdc++fs)
else()
        target_link_libraries(nano daedalum ${ITK_LIBRARIES})
endif()

