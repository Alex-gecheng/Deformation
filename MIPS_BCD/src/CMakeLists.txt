file(GLOB source
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

find_package(Eigen3 CONFIG REQUIRED)
add_executable(${PROJECT_NAME} ${source})

target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(${PROJECT_NAME} PRIVATE
    Eigen3::Eigen
)