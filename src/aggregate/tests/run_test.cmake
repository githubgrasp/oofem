# run_test.cmake — diff-based aggregate test driver
#
# Required arguments (pass via -D):
#   AGGREGATE   — path to aggregate.exec binary
#   TEST_DIR    — path to the test case directory
#
# Diffs `packing.dat` against `reference.dat` (always required). If a
# `reference.vtu` is present in the test directory, also diffs `packing.vtu`
# against it — this lets a single ctest entry cover both packing format and
# visualisation regression.

execute_process(
    COMMAND ${AGGREGATE} ${TEST_DIR}/control.in
    WORKING_DIRECTORY ${TEST_DIR}
    RESULT_VARIABLE run_result
)

if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "Aggregate exited with code ${run_result}")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files
        ${TEST_DIR}/packing.dat
        ${TEST_DIR}/reference.dat
    RESULT_VARIABLE diff_result
)

if(NOT diff_result EQUAL 0)
    message(FATAL_ERROR "Output 'packing.dat' differs from reference.dat in ${TEST_DIR}")
endif()

if(EXISTS ${TEST_DIR}/reference.vtu)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E compare_files
            ${TEST_DIR}/packing.vtu
            ${TEST_DIR}/reference.vtu
        RESULT_VARIABLE vtu_diff_result
    )
    if(NOT vtu_diff_result EQUAL 0)
        message(FATAL_ERROR "Output 'packing.vtu' differs from reference.vtu in ${TEST_DIR}")
    endif()
endif()
