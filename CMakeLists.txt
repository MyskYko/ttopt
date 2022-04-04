cmake_minimum_required(VERSION 3.8)
project(ttopt CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-g -O3 -Wall -DNDEBUG)

file(GLOB FILENAMES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(ttopt ${FILENAMES})
#target_include_directories(ttopt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
