cmake_minimum_required(VERSION 3.17)
project(optics-forward LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(OpenMP)

set(APP_SOURCES
  src/main.cpp
  src/tests.cpp
  src/utils.cpp
  src/surface.cpp
  src/surface_host.cpp
  src/poly_surface.cpp
  src/dep_definitions.cpp
  ext/tinyexr/miniz.c
)

add_executable(optics-forward ${APP_SOURCES})

target_include_directories(optics-forward PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/LiteMath
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/json/include
  ${CMAKE_CURRENT_SOURCE_DIR}/ext/tinyexr
  ${CMAKE_CURRENT_SOURCE_DIR}/ext
)

# The vendored dependencies contain warnings of their own; the tracing code is
# compiled with the project's normal warning level.
if(MSVC)
  set_source_files_properties(
    ext/tinyexr/miniz.c
    src/dep_definitions.cpp
    PROPERTIES COMPILE_OPTIONS "/w"
  )
else()
  set_source_files_properties(
    ext/tinyexr/miniz.c
    src/dep_definitions.cpp
    PROPERTIES COMPILE_OPTIONS "-w"
  )
endif()

if(OpenMP_CXX_FOUND)
  target_link_libraries(optics-forward PRIVATE OpenMP::OpenMP_CXX)
endif()

if(MSVC)
  target_compile_options(optics-forward PRIVATE /W3)
else()
  target_compile_options(optics-forward PRIVATE -Wall -Wextra -Wno-unused-parameter)
endif()

enable_testing()

add_test(
  NAME forward_tracing
  COMMAND optics-forward --output ${CMAKE_CURRENT_BINARY_DIR}/test-output
)

# Tests load 103.json and ref_data/* relative to the project root.
set_tests_properties(forward_tracing PROPERTIES
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
