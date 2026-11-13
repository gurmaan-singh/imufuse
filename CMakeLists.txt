cmake_minimum_required(VERSION 3.24)
project(imufuse VERSION 1.0.0 DESCRIPTION "IMU attitude estimation: complementary, Mahony and a multiplicative EKF" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(IMUFUSE_BUILD_TESTS "Build unit tests" ON)
option(IMUFUSE_WERROR "Treat warnings as errors" OFF)
option(IMUFUSE_SANITIZE "AddressSanitizer + UBSan (GCC/Clang)" OFF)

add_library(imufuse_warnings INTERFACE)
if(MSVC AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(imufuse_warnings INTERFACE /W4 /permissive- /utf-8 $<$<BOOL:${IMUFUSE_WERROR}>:/WX>)
else()
  target_compile_options(imufuse_warnings INTERFACE -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wold-style-cast
    $<$<BOOL:${IMUFUSE_WERROR}>:-Werror>)
endif()
if(IMUFUSE_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()

# header-only library
add_library(imufuse INTERFACE)
target_include_directories(imufuse INTERFACE include)

add_executable(imufuse_sim apps/imufuse_sim.cpp)
target_link_libraries(imufuse_sim PRIVATE imufuse imufuse_warnings)
target_compile_definitions(imufuse_sim PRIVATE $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>)

if(IMUFUSE_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
