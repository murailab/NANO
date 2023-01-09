message(STATUS "----------------Running CMake for Marching---------------------------")

target_sources(daedalum
        PRIVATE
            src/bwdist.cpp
            src/marching.cpp
            src/distances.cpp
        PUBLIC
            include/marching.h)

target_include_directories(daedalum PUBLIC include ${CMAKE_CURRENT_LIST_DIR})
if((${CMAKE_COMPILER_IS_GNUCXX}) AND (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS_EQUAL 8))
    target_link_libraries(daedalum PUBLIC  ${ITK_LIBRARIES} stdc++fs)
else()
    target_link_libraries(daedalum PUBLIC  ${ITK_LIBRARIES})
endif()
