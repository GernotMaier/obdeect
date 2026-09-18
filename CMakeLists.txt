cmake_minimum_required(VERSION 3.20)
project(obdeect VERSION 0.1.0 LANGUAGES CXX)

option(OBDEECT_BUILD_TESTING "Build obdeect tests" ON)

add_library(obdeect_core INTERFACE)
target_include_directories(obdeect_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include)
target_compile_features(obdeect_core INTERFACE cxx_std_20)

function(obdeect_enable_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
endfunction()

add_executable(obdeect_toy cpp/src/toy_main.cpp)
target_link_libraries(obdeect_toy PRIVATE obdeect_core)
obdeect_enable_warnings(obdeect_toy)

if(OBDEECT_BUILD_TESTING)
  enable_testing()
  add_executable(obdeect_core_tests cpp/tests/test_core.cpp)
  target_link_libraries(obdeect_core_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_core_tests)
  add_test(NAME obdeect_core_tests COMMAND obdeect_core_tests)

  add_executable(obdeect_literature_tests cpp/tests/test_literature_mst.cpp)
  target_link_libraries(obdeect_literature_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_literature_tests)
  add_test(NAME obdeect_literature_tests COMMAND obdeect_literature_tests)
endif()
