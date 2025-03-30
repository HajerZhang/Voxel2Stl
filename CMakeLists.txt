cmake_minimum_required(VERSION 3.20)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS OFF)

# message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
# message(STATUS "project source dir: ${CMAKE_SOURCE_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

project(
    Voxel2Stl
    VERSION 3.0
    LANGUAGES CXX C
)

add_subdirectory(src)

add_executable(
    main
    ${CMAKE_SOURCE_DIR}/src/main.cpp
)

target_link_libraries(
    main
    PRIVATE
    ${PROJECT_NAME}
)