find_package(Eigen3 REQUIRED NO_MODULE)

set(AGGREGATE_SOURCES
    box.C
    ellipsoid.C
    fibre.C
    intersection.C
    dlim.C
    gradingcurve.C
    placer.C
    vtu.C
   )

add_library(aggregate STATIC ${AGGREGATE_SOURCES})
target_link_libraries(aggregate PUBLIC liboofem Eigen3::Eigen)
target_include_directories(aggregate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(aggregate.exec aggregate.C)
target_link_libraries(aggregate.exec PRIVATE aggregate liboofem)

add_executable(intersection_test.exec intersection_test.C)
target_link_libraries(intersection_test.exec PRIVATE aggregate liboofem)

add_test(
    NAME test_aggregate_intersection_unit
    COMMAND $<TARGET_FILE:intersection_test.exec>
)

add_executable(grading_test.exec grading_test.C)
target_link_libraries(grading_test.exec PRIVATE aggregate liboofem)

add_test(
    NAME test_aggregate_grading_unit
    COMMAND $<TARGET_FILE:grading_test.exec>
)

add_executable(placer_test.exec placer_test.C)
target_link_libraries(placer_test.exec PRIVATE aggregate liboofem)

add_test(
    NAME test_aggregate_placer_unit
    COMMAND $<TARGET_FILE:placer_test.exec>
)

enable_testing()

# Aggregate tests — diff-based, one subdirectory per test case.
# Each subdirectory must contain control.in and reference.dat. The driver
# runs aggregate.exec on control.in and diffs the resulting packing against
# reference.dat.

set(AGGREGATE_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set(AGGREGATE_RUN_TEST ${AGGREGATE_TEST_DIR}/run_test.cmake)

file(GLOB _agg_test_entries RELATIVE ${AGGREGATE_TEST_DIR} ${AGGREGATE_TEST_DIR}/*)

foreach(_entry ${_agg_test_entries})
    set(_dir ${AGGREGATE_TEST_DIR}/${_entry})
    if(NOT IS_DIRECTORY ${_dir})
        continue()
    endif()
    if(NOT EXISTS ${_dir}/control.in OR NOT EXISTS ${_dir}/reference.dat)
        continue()
    endif()

    add_test(
        NAME test_aggregate_${_entry}
        COMMAND ${CMAKE_COMMAND}
            -DAGGREGATE=$<TARGET_FILE:aggregate.exec>
            -DTEST_DIR=${_dir}
            -P ${AGGREGATE_RUN_TEST}
    )
endforeach()
