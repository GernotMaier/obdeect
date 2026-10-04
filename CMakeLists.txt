cmake_minimum_required(VERSION 3.20)
project(obdeect VERSION 0.1.0 LANGUAGES CXX)

option(OBDEECT_BUILD_TESTING "Build obdeect tests" ON)
option(OBDEECT_BUILD_EVENTIO_INPUT "Build the CORSIKA EventIO photon reader" OFF)
if(OBDEECT_BUILD_EVENTIO_INPUT)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  find_package(Eventio REQUIRED)
endif()

if(DEFINED SKBUILD_PLATLIB_DIR)
  set(OBDEECT_BUILD_PYTHON_API_DEFAULT ON)
else()
  set(OBDEECT_BUILD_PYTHON_API_DEFAULT OFF)
endif()
option(OBDEECT_BUILD_PYTHON_API "Build the nanobind bulk tracing boundary"
       ${OBDEECT_BUILD_PYTHON_API_DEFAULT})

add_library(obdeect_core INTERFACE)
if(DEFINED SKBUILD_PLATLIB_DIR)
  set(OBDEECT_INSTALL_INCLUDE_DIR "obdeect/_include")
  set(OBDEECT_INSTALL_CMAKE_DIR "obdeect/cmake")
else()
  set(OBDEECT_INSTALL_INCLUDE_DIR "include")
  set(OBDEECT_INSTALL_CMAKE_DIR "lib/cmake/obdeect")
endif()
target_include_directories(obdeect_core INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/include>
  $<INSTALL_INTERFACE:${OBDEECT_INSTALL_INCLUDE_DIR}>)
target_compile_features(obdeect_core INTERFACE cxx_std_20)

function(obdeect_enable_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
endfunction()

if(OBDEECT_BUILD_EVENTIO_INPUT)
  add_library(obdeect_eventio cpp/src/eventio_photon_input.cpp)
  target_link_libraries(obdeect_eventio PUBLIC obdeect_core PRIVATE Eventio::Eventio)
  obdeect_enable_warnings(obdeect_eventio)
  add_executable(obdeect_eventio_summary cpp/src/eventio_summary_main.cpp)
  target_link_libraries(obdeect_eventio_summary PRIVATE obdeect_eventio)
  obdeect_enable_warnings(obdeect_eventio_summary)
endif()

add_executable(obdeect_demo_mst cpp/src/reference_main.cpp)
set_target_properties(obdeect_demo_mst PROPERTIES OUTPUT_NAME obdeect-demo-mst)
target_link_libraries(obdeect_demo_mst PRIVATE obdeect_core)
obdeect_enable_warnings(obdeect_demo_mst)

# Analytic prescriptions are a developer diagnostic, not CTAO telescope models.
add_executable(obdeect_analytic_optics cpp/src/ctao_main.cpp)
set_target_properties(obdeect_analytic_optics PROPERTIES OUTPUT_NAME obdeect-analytic-optics)
target_link_libraries(obdeect_analytic_optics PRIVATE obdeect_core)
obdeect_enable_warnings(obdeect_analytic_optics)

add_executable(obdeect-simtools-raytrace cpp/src/simtools_raytrace_main.cpp)
target_link_libraries(obdeect-simtools-raytrace PRIVATE obdeect_core)
obdeect_enable_warnings(obdeect-simtools-raytrace)

if(DEFINED SKBUILD_PLATLIB_DIR)
  # Wheels keep the native programs next to the Python launcher.  The launcher
  # resolves this directory through importlib, so users do not need a source
  # checkout or a separately installed runtime image.
  install(TARGETS obdeect_demo_mst obdeect_analytic_optics obdeect-simtools-raytrace
          RUNTIME DESTINATION obdeect/_native)
  install(DIRECTORY cpp/include/ DESTINATION obdeect/_include)
else()
  install(TARGETS obdeect_demo_mst obdeect_analytic_optics obdeect-simtools-raytrace RUNTIME DESTINATION bin)
  install(DIRECTORY cpp/include/ DESTINATION include)
endif()

if(OBDEECT_BUILD_EVENTIO_INPUT)
  if(DEFINED SKBUILD_PLATLIB_DIR)
    message(FATAL_ERROR "EventIO is a source-build capability; portable wheels disable it")
  endif()
  install(TARGETS obdeect_eventio_summary RUNTIME DESTINATION bin)
endif()

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

include(CMakePackageConfigHelpers)
install(TARGETS obdeect_core EXPORT obdeectTargets)
install(EXPORT obdeectTargets
        FILE obdeectTargets.cmake
        NAMESPACE obdeect::
        DESTINATION ${OBDEECT_INSTALL_CMAKE_DIR})
configure_package_config_file(
  cmake/obdeectConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/obdeectConfig.cmake
  INSTALL_DESTINATION ${OBDEECT_INSTALL_CMAKE_DIR})
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/obdeectConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/obdeectConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/obdeectConfigVersion.cmake
  DESTINATION ${OBDEECT_INSTALL_CMAKE_DIR})

if(OBDEECT_BUILD_TESTING)
  enable_testing()
  if(NOT DEFINED SKBUILD_PLATLIB_DIR)
    add_test(NAME obdeect_installed_consumer
      COMMAND ${CMAKE_COMMAND}
      -DOBDEECT_BUILD_DIR=${CMAKE_BINARY_DIR}
      -DOBDEECT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/tools/check_cmake_package.cmake)
  endif()
  # Keep test executables out of the user-facing binary directory.
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_TESTS "${CMAKE_BINARY_DIR}/tests")
  function(obdeect_add_test_target target source)
    add_executable(${target} ${source})
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_TESTS}")
    target_link_libraries(${target} PRIVATE obdeect_core)
    obdeect_enable_warnings(${target})
    # Scientific assertions must run in Release builds as well as Debug.
    if(MSVC)
      target_compile_options(${target} PRIVATE /UNDEBUG)
    else()
      target_compile_options(${target} PRIVATE -UNDEBUG)
    endif()
    add_test(NAME ${target} COMMAND ${target} ${ARGN})
  endfunction()
  obdeect_add_test_target(obdeect_detector_planes_tests cpp/tests/test_detector_planes.cpp)
  obdeect_add_test_target(obdeect_response_tests cpp/tests/test_optical_response.cpp)
  obdeect_add_test_target(obdeect_material_tests cpp/tests/test_materials.cpp)
  obdeect_add_test_target(obdeect_refractive_transport_tests cpp/tests/test_refractive_transport.cpp)
  obdeect_add_test_target(obdeect_core_tests cpp/tests/test_core.cpp)

  obdeect_add_test_target(obdeect_literature_tests cpp/tests/test_literature_mst.cpp)

  obdeect_add_test_target(obdeect_axisymmetric_optics_tests cpp/tests/test_axisymmetric_optics.cpp)

  obdeect_add_test_target(obdeect_source_sampling_tests cpp/tests/test_source_sampling.cpp)
  obdeect_add_test_target(obdeect_source_tests cpp/tests/test_sources.cpp)

  obdeect_add_test_target(obdeect_photon_input_tests cpp/tests/test_photon_input.cpp)

  obdeect_add_test_target(obdeect_atmosphere_tests cpp/tests/test_atmosphere.cpp)

  if(OBDEECT_BUILD_EVENTIO_INPUT)
    add_executable(obdeect_eventio_input_tests cpp/tests/test_eventio_photon_input.cpp)
    set_target_properties(obdeect_eventio_input_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY_TESTS}")
    target_link_libraries(obdeect_eventio_input_tests PRIVATE obdeect_eventio Eventio::Eventio)
    obdeect_enable_warnings(obdeect_eventio_input_tests)
    if(MSVC)
      target_compile_options(obdeect_eventio_input_tests PRIVATE /UNDEBUG)
    else()
      target_compile_options(obdeect_eventio_input_tests PRIVATE -UNDEBUG)
    endif()
    add_test(NAME obdeect_eventio_input_tests COMMAND obdeect_eventio_input_tests)
  endif()

  obdeect_add_test_target(obdeect_native_replay_tests cpp/tests/test_native_replay.cpp
                         $<TARGET_FILE:obdeect-simtools-raytrace>)
  obdeect_add_test_target(obdeect_trace_tests cpp/tests/test_trace.cpp)

  obdeect_add_test_target(obdeect_contract_tests cpp/tests/test_contracts.cpp)

  obdeect_add_test_target(obdeect_ctao_model_tests cpp/tests/test_ctao_models.cpp)

  obdeect_add_test_target(obdeect_facet_tests cpp/tests/test_facets.cpp)

  obdeect_add_test_target(obdeect_segmented_optical_model_tests cpp/tests/test_segmented_optical_model.cpp)

  obdeect_add_test_target(obdeect_optical_model_file_tests cpp/tests/test_optical_model_file.cpp)

  obdeect_add_test_target(obdeect_optical_model_tests cpp/tests/test_optical_model.cpp)

endif()
