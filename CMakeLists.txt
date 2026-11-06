cmake_minimum_required(VERSION 3.20)
project(constrained_xpbd LANGUAGES CXX)

find_path(EIGEN_INCLUDE_DIR Eigen/Core HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external/eigen")
if(NOT EIGEN_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen missing. Run git submodule update --init external/eigen")
endif()

add_library(xpbd INTERFACE)
target_compile_features(xpbd INTERFACE cxx_std_17)
target_include_directories(xpbd INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(xpbd SYSTEM INTERFACE "${EIGEN_INCLUDE_DIR}")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(xpbd INTERFACE -Wall -Wextra -Wpedantic)
endif()

add_executable(constrained-xpbd main.cpp)
target_link_libraries(constrained-xpbd PRIVATE xpbd)

include(CTest)
if(BUILD_TESTING)
    add_executable(xpbd-baseline-tests tests.cpp)
    target_link_libraries(xpbd-baseline-tests PRIVATE xpbd)
    add_test(NAME baseline-and-passive COMMAND xpbd-baseline-tests)
    add_executable(xpbd-active-tests active_tests.cpp)
    target_link_libraries(xpbd-active-tests PRIVATE xpbd)
    add_test(NAME active-reservoirs COMMAND xpbd-active-tests)
    set_tests_properties(active-reservoirs PROPERTIES DEPENDS baseline-and-passive)
    add_executable(xpbd-sequential-tests sequential_tests.cpp)
    target_link_libraries(xpbd-sequential-tests PRIVATE xpbd)
    add_test(NAME sequential-work-feedback COMMAND xpbd-sequential-tests)
    set_tests_properties(sequential-work-feedback PROPERTIES DEPENDS baseline-and-passive)
endif()

# Research probes: nonlinear work balance, curvature-dropped linearization, dual GS.
option(XPBD_BUILD_WORK_STUDIES "Build the spring/pendulum work-balance studies" ON)
if(XPBD_BUILD_WORK_STUDIES)
    foreach(study IN ITEMS kinetic_work midpoint_gs spring pendulum coupled_momentum)
        string(REPLACE "_" "-" study_name "${study}")
        add_executable(xpbd-${study_name}-study experiments/work_balance/${study}_study.cpp)
        target_link_libraries(xpbd-${study_name}-study PRIVATE xpbd)
        if(BUILD_TESTING)
            file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/work_balance/${study}")
            add_test(NAME work-${study_name} COMMAND xpbd-${study_name}-study)
            set_tests_properties(work-${study_name} PROPERTIES
                WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/work_balance/${study}"
                DEPENDS baseline-and-passive)
        endif()
    endforeach()
endif()

# Direct energy conservation: joint multiplier blocks versus sequential updates.
option(XPBD_BUILD_ENERGY_STUDIES "Build the direct-energy coupling validation" ON)
if(XPBD_BUILD_ENERGY_STUDIES)
    foreach(study IN ITEMS block spring pendulum)
        add_executable(xpbd-energy-${study}-study experiments/energy_coupling/${study}_study.cpp)
        target_link_libraries(xpbd-energy-${study}-study PRIVATE xpbd)
        if(BUILD_TESTING)
            file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/energy_coupling/${study}")
            add_test(NAME energy-coupling-${study} COMMAND xpbd-energy-${study}-study)
            set_tests_properties(energy-coupling-${study} PROPERTIES
                WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/energy_coupling/${study}"
                DEPENDS baseline-and-passive)
        endif()
    endforeach()
endif()

# Compare paired and sequential local work feedback with separate compliances.
if(XPBD_BUILD_WORK_STUDIES)
    add_executable(xpbd-ordering-study experiments/sequential_work/ordering_study.cpp)
    target_link_libraries(xpbd-ordering-study PRIVATE xpbd)
endif()

# Six global-energy schedules, with compliant energy and retained-r_x controls.
if(XPBD_BUILD_ENERGY_STUDIES)
    add_executable(xpbd-energy-ordering-study experiments/energy_ordering/study.cpp)
    target_link_libraries(xpbd-energy-ordering-study PRIVATE xpbd)
    if(BUILD_TESTING)
        add_executable(xpbd-energy-ordering-tests experiments/energy_ordering/tests.cpp)
        target_link_libraries(xpbd-energy-ordering-tests PRIVATE xpbd)
        add_test(NAME energy-ordering COMMAND xpbd-energy-ordering-tests)
        set_tests_properties(energy-ordering PROPERTIES DEPENDS baseline-and-passive)
    endif()
endif()

# Compare GS orderings of one frozen material/energy linear system.
if(XPBD_BUILD_ENERGY_STUDIES)
    add_executable(xpbd-frozen-schedules-study experiments/energy_ordering/frozen_study.cpp)
    target_link_libraries(xpbd-frozen-schedules-study PRIVATE xpbd)
    if(BUILD_TESTING)
        add_executable(xpbd-frozen-schedules-tests experiments/energy_ordering/frozen_tests.cpp)
        target_link_libraries(xpbd-frozen-schedules-tests PRIVATE xpbd)
        add_test(NAME frozen-energy-schedules COMMAND xpbd-frozen-schedules-tests)
        set_tests_properties(frozen-energy-schedules PROPERTIES DEPENDS baseline-and-passive)
    endif()
endif()

# Spring diagnosis: direct-energy controls and a midpoint XPBD comparison.
if(XPBD_BUILD_ENERGY_STUDIES)
    add_executable(xpbd-spring-midpoint-study experiments/spring_repair/midpoint_study.cpp)
    target_link_libraries(xpbd-spring-midpoint-study PRIVATE xpbd)
    add_executable(xpbd-spring-work-midpoint-study experiments/spring_repair/work_midpoint_study.cpp)
    target_link_libraries(xpbd-spring-work-midpoint-study PRIVATE xpbd)
    add_executable(xpbd-spring-controls-study experiments/spring_repair/controls_study.cpp)
    target_link_libraries(xpbd-spring-controls-study PRIVATE xpbd)
    target_include_directories(xpbd-spring-controls-study PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/experiments/energy_ordering")
    if(BUILD_TESTING)
        add_test(NAME spring-midpoint-validation COMMAND xpbd-spring-midpoint-study
            "${CMAKE_BINARY_DIR}/spring_midpoint_validation")
        set_tests_properties(spring-midpoint-validation PROPERTIES DEPENDS baseline-and-passive)
        add_test(NAME spring-work-midpoint-validation COMMAND xpbd-spring-work-midpoint-study
            "${CMAKE_BINARY_DIR}/spring_work_midpoint_validation")
        set_tests_properties(spring-work-midpoint-validation PROPERTIES DEPENDS baseline-and-passive)
    endif()
endif()

# Original endpoint-energy equations: oversized update rejection and old gradients.
if(XPBD_BUILD_ENERGY_STUDIES)
    add_executable(xpbd-gradient-policies-study experiments/gradient_policies/policy_study.cpp)
    target_link_libraries(xpbd-gradient-policies-study PRIVATE xpbd)
    target_include_directories(xpbd-gradient-policies-study PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/experiments/energy_ordering")
    if(BUILD_TESTING)
        add_executable(xpbd-gradient-policies-tests experiments/gradient_policies/policy_tests.cpp)
        target_link_libraries(xpbd-gradient-policies-tests PRIVATE xpbd)
        target_include_directories(xpbd-gradient-policies-tests PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/experiments/energy_ordering")
        add_test(NAME gradient-policies COMMAND xpbd-gradient-policies-tests)
        set_tests_properties(gradient-policies PROPERTIES DEPENDS baseline-and-passive)
    endif()
endif()

# Gravity plus a compliant spring: validate EL before extending to a chain.
if(XPBD_BUILD_ENERGY_STUDIES)
    add_executable(xpbd-elastic-pendulum-study experiments/elastic_pendulum/study.cpp)
    target_link_libraries(xpbd-elastic-pendulum-study PRIVATE xpbd)
    target_include_directories(xpbd-elastic-pendulum-study PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/experiments/gradient_policies"
        "${CMAKE_CURRENT_SOURCE_DIR}/experiments/energy_ordering")
    if(BUILD_TESTING)
        add_executable(xpbd-elastic-pendulum-tests experiments/elastic_pendulum/tests.cpp)
        target_link_libraries(xpbd-elastic-pendulum-tests PRIVATE xpbd)
        target_include_directories(xpbd-elastic-pendulum-tests PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/experiments/gradient_policies"
            "${CMAKE_CURRENT_SOURCE_DIR}/experiments/energy_ordering")
        add_test(NAME elastic-pendulum COMMAND xpbd-elastic-pendulum-tests)
        set_tests_properties(elastic-pendulum PROPERTIES DEPENDS baseline-and-passive)
    endif()
endif()

# Reject energy increments that oppose descent of the primal implicit objective.
if(XPBD_BUILD_ENERGY_STUDIES AND BUILD_TESTING)
    add_executable(xpbd-directional-clip-tests experiments/elastic_pendulum/direction_tests.cpp)
    target_link_libraries(xpbd-directional-clip-tests PRIVATE xpbd)
    target_include_directories(xpbd-directional-clip-tests PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/experiments/gradient_policies"
        "${CMAKE_CURRENT_SOURCE_DIR}/experiments/energy_ordering")
    add_test(NAME directional-energy-gate COMMAND xpbd-directional-clip-tests)
    set_tests_properties(directional-energy-gate PROPERTIES DEPENDS baseline-and-passive)
endif()

# Small steps: one frozen EL sweep with damping, capping, or compliant energy.
if(XPBD_BUILD_ENERGY_STUDIES)
    add_executable(xpbd-small-steps-study experiments/small_steps/study.cpp)
    target_link_libraries(xpbd-small-steps-study PRIVATE xpbd)
    target_include_directories(xpbd-small-steps-study PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/experiments/gradient_policies"
        "${CMAKE_CURRENT_SOURCE_DIR}/experiments/energy_ordering")
    if(BUILD_TESTING)
        add_executable(xpbd-small-steps-tests experiments/small_steps/tests.cpp)
        target_link_libraries(xpbd-small-steps-tests PRIVATE xpbd)
        target_include_directories(xpbd-small-steps-tests PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/experiments/gradient_policies"
            "${CMAKE_CURRENT_SOURCE_DIR}/experiments/energy_ordering")
        add_test(NAME small-steps-el-controls COMMAND xpbd-small-steps-tests)
        set_tests_properties(small-steps-el-controls PROPERTIES DEPENDS baseline-and-passive)
    endif()
endif()
