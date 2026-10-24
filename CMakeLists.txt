cmake_minimum_required(VERSION 3.20)
project(obdeect VERSION 0.2.0 LANGUAGES CXX)

option(OBDEECT_BUILD_TESTING "Build the C++ tests" ON)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(OBDEECT_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(OBDEECT_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(OBDEECT_BUILD_BENCH "Build the benchmark and replay tools (not run by ctest)" OFF)
option(OBDEECT_BUILD_PYTHON_API "Build the Python extension" ${OBDEECT_BUILD_PYTHON_API_DEFAULT})

# Header-only, standard-library-only core.
add_library(obdeect_core INTERFACE)
target_include_directories(obdeect_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_compile_features(obdeect_core INTERFACE cxx_std_20)
# Identical floating-point results on every platform: no implicit fused
# multiply-add contraction.
if(NOT MSVC)
  target_compile_options(obdeect_core INTERFACE -ffp-contract=off)
endif()

function(obdeect_enable_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
endfunction()

if(OBDEECT_BUILD_PYTHON_API)
  find_package(Python 3.14 COMPONENTS Interpreter Development.Module REQUIRED)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
                  OUTPUT_VARIABLE nanobind_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
                  COMMAND_ERROR_IS_FATAL ANY)
  find_package(nanobind CONFIG REQUIRED)
  nanobind_add_module(_core NB_STATIC cpp/bindings/core.cpp)
  target_link_libraries(_core PRIVATE obdeect_core)
  obdeect_enable_warnings(_core)
  install(TARGETS _core LIBRARY DESTINATION obdeect)
endif()

if(OBDEECT_BUILD_TESTING)
  enable_testing()
  foreach(name kernels trace model_json)
    add_executable(test_${name} cpp/tests/test_${name}.cpp)
    target_link_libraries(test_${name} PRIVATE obdeect_core)
    obdeect_enable_warnings(test_${name})
    add_test(NAME ${name} COMMAND test_${name})
  endforeach()
endif()

if(OBDEECT_BUILD_BENCH)
  foreach(name bench_trace replay)
    add_executable(${name} cpp/bench/${name}.cpp)
    target_link_libraries(${name} PRIVATE obdeect_core)
    obdeect_enable_warnings(${name})
  endforeach()
endif()
