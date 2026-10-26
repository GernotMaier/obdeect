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
# Identical floating-point results for any thread count and block size on one
# platform: no implicit fused multiply-add contraction. Transcendental functions
# come from the platform's math library, so results may differ between
# platforms.
if(NOT MSVC)
  target_compile_options(obdeect_core INTERFACE -ffp-contract=off)
endif()

# Changes results at the rounding level. Use only on hardware with fused
# multiply-add; without it std::fma is emulated and very slow. Results differ
# from builds without it. Default builds, wheels and CI stay without FMA.
option(OBDEECT_USE_FMA "Use std::fma in the aspheric polynomial evaluation" OFF)
if(OBDEECT_USE_FMA)
  target_compile_definitions(obdeect_core INTERFACE OBDEECT_USE_FMA=1)
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
  foreach(name kernels trace model_json bounds)
    add_executable(test_${name} cpp/tests/test_${name}.cpp)
    target_link_libraries(test_${name} PRIVATE obdeect_core)
    obdeect_enable_warnings(test_${name})
    add_test(NAME ${name} COMMAND test_${name})
  endforeach()
  # Newton-step counts need the benchmark counters compiled in.
  add_executable(test_asphere_steps cpp/tests/test_asphere_steps.cpp)
  target_link_libraries(test_asphere_steps PRIVATE obdeect_core)
  target_compile_definitions(test_asphere_steps PRIVATE OBDEECT_BENCH_COUNTERS)
  obdeect_enable_warnings(test_asphere_steps)
  add_test(NAME asphere_steps COMMAND test_asphere_steps)
  # The replay tool's limit comparison (header in cpp/bench, standard library only).
  add_executable(test_replay_compare cpp/tests/test_replay_compare.cpp)
  target_link_libraries(test_replay_compare PRIVATE obdeect_core)
  target_include_directories(test_replay_compare PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/bench)
  obdeect_enable_warnings(test_replay_compare)
  add_test(NAME replay_compare COMMAND test_replay_compare)
endif()

if(OBDEECT_BUILD_BENCH)
  foreach(name bench_trace replay)
    add_executable(${name} cpp/bench/${name}.cpp)
    target_link_libraries(${name} PRIVATE obdeect_core)
    obdeect_enable_warnings(${name})
  endforeach()
endif()
