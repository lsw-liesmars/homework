cmake_minimum_required(VERSION 3.10)
project(PnP_Solver)

# 设置C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2")

# 主程序（包含所有代码）
add_executable(pnp_solver pnp_solver.cpp)

# 测试程序（注意：由于test_pnp.cpp包含了pnp_solver.cpp，直接编译即可）
# 如果编译器报错重复定义，可以将pnp_solver.cpp改为头文件形式
add_executable(test_pnp test_pnp.cpp)

# 链接数学库
target_link_libraries(pnp_solver m)
target_link_libraries(test_pnp m)

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS}")
