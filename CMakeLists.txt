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

add_executable(obdeect_ctao cpp/src/ctao_main.cpp)
target_link_libraries(obdeect_ctao PRIVATE obdeect_core)
obdeect_enable_warnings(obdeect_ctao)

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

  add_executable(obdeect_axisymmetric_optics_tests cpp/tests/test_axisymmetric_optics.cpp)
  target_link_libraries(obdeect_axisymmetric_optics_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_axisymmetric_optics_tests)
  add_test(NAME obdeect_axisymmetric_optics_tests COMMAND obdeect_axisymmetric_optics_tests)

  add_executable(obdeect_source_tests cpp/tests/test_sources.cpp)
  target_link_libraries(obdeect_source_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_source_tests)
  add_test(NAME obdeect_source_tests COMMAND obdeect_source_tests)

  add_executable(obdeect_trace_tests cpp/tests/test_trace.cpp)
  target_link_libraries(obdeect_trace_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_trace_tests)
  add_test(NAME obdeect_trace_tests COMMAND obdeect_trace_tests)

  add_executable(obdeect_contract_tests cpp/tests/test_contracts.cpp)
  target_link_libraries(obdeect_contract_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_contract_tests)
  add_test(NAME obdeect_contract_tests COMMAND obdeect_contract_tests)

  add_executable(obdeect_ctao_model_tests cpp/tests/test_ctao_models.cpp)
  target_link_libraries(obdeect_ctao_model_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_ctao_model_tests)
  add_test(NAME obdeect_ctao_model_tests COMMAND obdeect_ctao_model_tests)

  add_executable(obdeect_facet_tests cpp/tests/test_facets.cpp)
  target_link_libraries(obdeect_facet_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_facet_tests)
  add_test(NAME obdeect_facet_tests COMMAND obdeect_facet_tests)
endif()
