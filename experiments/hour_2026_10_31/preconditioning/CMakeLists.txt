cmake_minimum_required(VERSION 3.20)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(xpbd_hour_preconditioning LANGUAGES CXX)
  enable_testing()
endif()

# The experiment owns its model helpers; Eigen is the only checkout dependency.
set(XPBD_ROOT "$ENV{HOME}/Projects/constrained-xpbd" CACHE PATH
    "Constrained-XPBD checkout containing external/eigen")
find_path(XPBD_EIGEN_INCLUDE_DIR Eigen/Dense
          HINTS "${XPBD_ROOT}/external/eigen" /usr/include/eigen3 REQUIRED)

foreach(driver benchmark physical ablation curvature_graph forcing cache refinement timing_followup)
  add_executable(xpbd_hour_preconditioning_${driver} ${driver}.cpp)
  target_compile_features(xpbd_hour_preconditioning_${driver} PRIVATE cxx_std_17)
  target_include_directories(xpbd_hour_preconditioning_${driver} PRIVATE
                             "${XPBD_EIGEN_INCLUDE_DIR}")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(xpbd_hour_preconditioning_${driver} PRIVATE
                           $<$<CONFIG:Release>:-O2>)
  endif()
endforeach()

# Write validation CSVs into the build tree, never over the recorded experiments.
foreach(check identity curvature_graph cached_factors)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/validation/${check}")
endforeach()
add_test(NAME preconditioning.identity
         COMMAND xpbd_hour_preconditioning_benchmark --tests)
set_tests_properties(preconditioning.identity PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/validation/identity"
  TIMEOUT 120)
add_test(NAME preconditioning.curvature_graph
         COMMAND xpbd_hour_preconditioning_curvature_graph)
set_tests_properties(preconditioning.curvature_graph PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/validation/curvature_graph"
  TIMEOUT 120)
add_test(NAME preconditioning.cached_factors
         COMMAND xpbd_hour_preconditioning_cache)
set_tests_properties(preconditioning.cached_factors PROPERTIES
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/validation/cached_factors"
  TIMEOUT 120)
