cmake_minimum_required(VERSION 3.16)
project(RBF LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# 查找 Eigen
find_package(Eigen3 3.3 REQUIRED)
message(STATUS "Eigen3 include dir: ${EIGEN3_INCLUDE_DIR}")

# 源码
file(GLOB_RECURSE SRC_FILES                                  
    ${PROJECT_SOURCE_DIR}/src/*.cpp
)

add_executable(${PROJECT_NAME} ${SRC_FILES})

# 工程头文件
target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${EIGEN3_INCLUDE_DIR}  # ✅ 添加 Eigen 头文件路径
)

# 链接库（Eigen 是头文件库，可选）
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
