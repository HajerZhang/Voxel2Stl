set(SRC_LIST
    voxel2stl.cpp
)

find_package(VTK REQUIRED)
message(STATUS "VTK_LIBRARIES: ${VTK_LIBRARIES}")

add_library(
    ${PROJECT_NAME} 
    ${SRC_LIST}
)

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    # ${VTK_USE_FILE} // VTK_USE_FILE is not used starting with VTK 8.90, only link VTK_LIBRARIES
)

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
    ${VTK_LIBRARIES}
)

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
)