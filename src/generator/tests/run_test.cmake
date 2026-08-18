# run_test.cmake — diff-based generator test driver
#
# Required arguments (pass via -D):
#   GENERATOR   — path to generator.exec binary
#   TEST_DIR    — path to the test case directory (contains mesh.in + reference.nodes)

execute_process(
    COMMAND ${GENERATOR} ${TEST_DIR}/mesh.in
    WORKING_DIRECTORY ${TEST_DIR}
    RESULT_VARIABLE run_result
)

if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "Generator exited with code ${run_result}")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files
        ${TEST_DIR}/nodes.dat
        ${TEST_DIR}/reference.nodes
    RESULT_VARIABLE diff_result
)

if(NOT diff_result EQUAL 0)
    message(FATAL_ERROR "Output 'nodes.dat' differs from reference.nodes in ${TEST_DIR}")
endif()
