cmake_minimum_required(VERSION 3.16)
project(poly2d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(poly2d src/main.cpp)
target_include_directories(poly2d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(MSVC)
  target_compile_options(poly2d PRIVATE /O2)
else()
  target_compile_options(poly2d PRIVATE -O2 -Wall -Wextra)
endif()
