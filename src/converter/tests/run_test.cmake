# run_test.cmake — diff-based converter test driver
#
# Required arguments (pass via -D):
#   CONVERTER   — path to converter.exec binary
#   TEST_DIR    — path to the test case directory
#   MESH_ARGS   — semicolon-separated list of mesh file paths passed after control.in
#   OUTPUT      — name of the output file produced by the converter

execute_process(
    COMMAND ${CONVERTER} ${TEST_DIR}/control.in ${MESH_ARGS}
    WORKING_DIRECTORY ${TEST_DIR}
    RESULT_VARIABLE run_result
)

if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "Converter exited with code ${run_result}")
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files
        ${TEST_DIR}/${OUTPUT}
        ${TEST_DIR}/reference.in
    RESULT_VARIABLE diff_result
)

if(NOT diff_result EQUAL 0)
    message(FATAL_ERROR "Output '${OUTPUT}' differs from reference.in in ${TEST_DIR}")
endif()
