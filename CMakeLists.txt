cmake_minimum_required(VERSION 3.20)
project(obdeect VERSION 0.1.0 LANGUAGES CXX)

option(OBDEECT_BUILD_TESTING "Build obdeect tests" ON)
option(OBDEECT_BUILD_EVENTIO_INPUT "Build the CORSIKA EventIO photon reader" OFF)
if(OBDEECT_BUILD_EVENTIO_INPUT)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  find_package(Eventio REQUIRED)
endif()

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

add_executable(obdeect_toy cpp/src/toy_main.cpp)
target_link_libraries(obdeect_toy PRIVATE obdeect_core)
obdeect_enable_warnings(obdeect_toy)

add_executable(obdeect_ctao cpp/src/ctao_main.cpp)
target_link_libraries(obdeect_ctao PRIVATE obdeect_core)
obdeect_enable_warnings(obdeect_ctao)

add_executable(obdeect-simtools-raytrace cpp/src/simtools_raytrace_main.cpp)
target_link_libraries(obdeect-simtools-raytrace PRIVATE obdeect_core)
obdeect_enable_warnings(obdeect-simtools-raytrace)

if(DEFINED SKBUILD_PLATLIB_DIR)
  # Wheels keep the native programs next to the Python launcher.  The launcher
  # resolves this directory through importlib, so users do not need a source
  # checkout or a separately installed runtime image.
  install(TARGETS obdeect_toy obdeect_ctao obdeect-simtools-raytrace
          RUNTIME DESTINATION obdeect/_native)
  install(DIRECTORY cpp/include/ DESTINATION obdeect/_include)
else()
  install(TARGETS obdeect_toy obdeect_ctao obdeect-simtools-raytrace RUNTIME DESTINATION bin)
  install(DIRECTORY cpp/include/ DESTINATION include)
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
  add_executable(obdeect_material_tests cpp/tests/test_materials.cpp)
  target_link_libraries(obdeect_material_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_material_tests)
  add_test(NAME obdeect_material_tests COMMAND obdeect_material_tests)
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

  add_executable(obdeect_photon_input_tests cpp/tests/test_photon_input.cpp)
  target_link_libraries(obdeect_photon_input_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_photon_input_tests)
  add_test(NAME obdeect_photon_input_tests COMMAND obdeect_photon_input_tests)

  add_executable(obdeect_atmosphere_tests cpp/tests/test_atmosphere.cpp)
  target_link_libraries(obdeect_atmosphere_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_atmosphere_tests)
  add_test(NAME obdeect_atmosphere_tests COMMAND obdeect_atmosphere_tests)

  if(OBDEECT_BUILD_EVENTIO_INPUT)
    add_executable(obdeect_eventio_input_tests cpp/tests/test_eventio_photon_input.cpp)
    target_link_libraries(obdeect_eventio_input_tests PRIVATE obdeect_eventio Eventio::Eventio)
    obdeect_enable_warnings(obdeect_eventio_input_tests)
    add_test(NAME obdeect_eventio_input_tests COMMAND obdeect_eventio_input_tests)
  endif()

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

  add_executable(obdeect_segmented_scene_tests cpp/tests/test_segmented_scene.cpp)
  target_link_libraries(obdeect_segmented_scene_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_segmented_scene_tests)
  add_test(NAME obdeect_segmented_scene_tests COMMAND obdeect_segmented_scene_tests)

  add_executable(obdeect_optical_scene_tests cpp/tests/test_optical_scene.cpp)
  target_link_libraries(obdeect_optical_scene_tests PRIVATE obdeect_core)
  obdeect_enable_warnings(obdeect_optical_scene_tests)
  add_test(NAME obdeect_optical_scene_tests COMMAND obdeect_optical_scene_tests)
endif()
