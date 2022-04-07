cmake_minimum_required(VERSION 3.8)
project(ttopt CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-g -O3 -Wall -DNDEBUG)

add_subdirectory(cudd)
target_include_directories(cudd INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/cudd/cudd;${CMAKE_CURRENT_SOURCE_DIR}/cudd/cplusplus")

file(GLOB FILENAMES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(ttopt ${FILENAMES})
target_link_libraries(ttopt cudd)
#target_include_directories(ttopt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
