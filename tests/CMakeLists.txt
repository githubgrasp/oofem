# tests/regression/CMakeLists.txt

# Ensure CTest is enabled
include(CTest)

# Find Python if some regression tests use Python scripts
find_package(Python3 COMPONENTS Interpreter)

if(USE_OOFEM_EXE)
set (oofem_REGRESSION_TEST_DIR "${oofem_SOURCE_DIR}/tests/regression")
set (oofem_BENCHMARK_DIR "${oofem_REGRESSION_TEST_DIR}/benchmark")
set (oofem_cmd $<TARGET_FILE:oofem>)
#file(TO_CMAKE_PATH "${oofem_TEST_DIR}/valgrind.supp" MEMORYCHECK_SUPPRESSIONS_FILE)
#set (CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${oofem_TEST_DIR}/tests/valgrind.supp)

if (USE_MPI_PARALLEL)
    set (test_rel_dir "tests/regression/partests")
    set (test_dir "${oofem_SOURCE_DIR}/${test_rel_dir}")
    file (GLOB parallel_tests RELATIVE "${test_dir}" "${test_dir}/*/")
    foreach (case ${parallel_tests})
        file (GLOB files "${test_dir}/${case}/${case}.oofem.in.*")
        list (LENGTH files num_files)
        add_test (NAME "${test_rel_dir}/${case}" WORKING_DIRECTORY ${test_dir}/${case} COMMAND "mpirun" "--oversubscribe" "-np" ${num_files} ${oofem_cmd} "-p" "-f" ${case}.oofem.in "-ksp_type" "cg")
    endforeach (case)
    # partest_brazil_2d_nl7 requires more time on testing server than default 1500
    set_tests_properties(${test_rel_dir}/brazil_2d_nl7 PROPERTIES TIMEOUT 2500)
endif ()

set (oofem_cmd $<TARGET_FILE:oofem> -Werror)
function(oofem_regression_ctest dir label)
    # discover tests
    set (test_rel_dir "tests/regression")
    file(GLOB gg RELATIVE "${oofem_SOURCE_DIR}/${test_rel_dir}/${dir}" "${oofem_SOURCE_DIR}/${test_rel_dir}/${dir}/*.in")
    foreach(case ${gg})
        add_test(NAME "${test_rel_dir}/${dir}/${case}" WORKING_DIRECTORY "${oofem_SOURCE_DIR}/${test_rel_dir}/${dir}" COMMAND ${oofem_cmd} "-f" ${case})
        set_tests_properties("${test_rel_dir}/${dir}/${case}" PROPERTIES LABELS ${label})
    endforeach(case)

    # discover benchmarks
    set (benchmark_rel_dir "tests/regression/benchmark")
    file(GLOB gg RELATIVE "${oofem_SOURCE_DIR}/${benchmark_rel_dir}/${dir}" "${oofem_BENCHMARK_DIR}/${dir}/*.in")
    foreach(case ${gg})
        add_test(NAME "${benchmark_rel_dir}/${dir}/${case}" WORKING_DIRECTORY "${oofem_SOURCE_DIR}/${benchmark_rel_dir}/${dir}" COMMAND ${oofem_cmd} "-f" ${case})
        set_tests_properties("${benchmark_rel_dir}/${dir}/${case}" PROPERTIES LABELS "${label};benchmark")
    endforeach(case)
    

    # repeat inside xml/ for both tests and benchmarks, if XML is enabled
    if(USE_XML)
        #file(GLOB gg RELATIVE "${oofem_TEST_DIR}/xml/${dir}" "${oofem_TEST_DIR}/xml/${dir}/*.xml")
        set (test_rel_dir "tests/regression")
        file(GLOB gg RELATIVE "${oofem_SOURCE_DIR}/${test_rel_dir}/${dir}" "${oofem_REGRESSION_TEST_DIR}/${dir}/*.xml")
        foreach(case ${gg})
            # skip subproblems; they could be called something else than *.xml but then DataReader would not detect XML correctly
            if(${case} MATCHES "\\.(fragment|frag)\\.xml$")
                continue()
            endif()
            add_test(NAME "${test_rel_dir}/${dir}/${case}" WORKING_DIRECTORY "${oofem_SOURCE_DIR}/${test_rel_dir}/${dir}" COMMAND ${oofem_cmd} "-f" ${case})
            set_tests_properties("${test_rel_dir}/${dir}/${case}" PROPERTIES LABELS ${label})
        endforeach(case)
    endif()
endfunction()

function(add_oofem_benchmarks dir wildcard)
endfunction()


# Sequential test cases
if (USE_SM)
    oofem_regression_ctest(sm "sm")
    oofem_regression_ctest(sm/NAFEMS "sm")


    if (NOT WIN32)
        message(STATUS "Identifying workflow tests for SM module 2")

        set (test_rel_dir "tests/regression/sm")
        file (GLOB sm_tests RELATIVE "${oofem_SOURCE_DIR}/${test_rel_dir}" "${oofem_SOURCE_DIR}/${test_rel_dir}/*.sh")
        foreach (case ${sm_tests})
            add_test (NAME "${test_rel_dir}/${case}" WORKING_DIRECTORY ${oofem_SOURCE_DIR}/${test_rel_dir} COMMAND bash ${case} ${oofem_cmd})
            set_tests_properties("${test_rel_dir}/${case}" PROPERTIES LABELS "workflow")
        endforeach (case)
    endif()
endif ()


if (USE_FM)
    oofem_regression_ctest(fm "fm")
endif()

if (USE_TM)
    oofem_regression_ctest(tm "tm")
endif()

if (USE_AM)
    oofem_regression_ctest(am "am")
endif()

if (USE_TM AND USE_SM)
    oofem_regression_ctest(tmsm "tmsm")
endif()

if (USE_TM AND USE_FM)
    oofem_regression_ctest(tmfm "tmfm")
endif()

if (USE_TM AND USE_CEMHYD)
    oofem_regression_ctest(tmcemhyd "tmcemhyd")
endif()

if (USE_FM AND USE_SM AND USE_PFEM)
    oofem_regression_ctest(fmpfem "fmpfem")
endif()

if (USE_SM AND USE_MFRONT)
    oofem_regression_ctest(smmfront "smmfront")
endif()

if (USE_SM AND USE_PYTHON_EXTENSION)
    oofem_regression_ctest(sm/python "sm/python")
endif()

if (USE_MPM)
    oofem_regression_ctest(mpm "mpm")
endif()

if (USE_TM AND USE_SM AND USE_DSS AND USE_IML)
    oofem_regression_ctest(tmsmdssiml "tmsmdssiml")
endif()


# add unit tests from tests/unit

add_subdirectory(unit)

# Custom targets for test and benchmarks (which can conveniently be accessed from IDEs)

#ProcessorCount(N)
if (N EQUAL 0)
  set (N 1)
endif ()
if (USE_OPENMP_PARALLEL) # OPENMP will use all processors by default
  set (N 1)
endif ()

if (MSVC)
    add_custom_target (tests      COMMAND ${CMAKE_CTEST_COMMAND} "-j${N}" "-R" "^test" "-C" "Debug")
    add_custom_target (partests   COMMAND ${CMAKE_CTEST_COMMAND} "-R" "^partest" "-C" "Debug")
    add_custom_target (benchmarks COMMAND ${CMAKE_CTEST_COMMAND} "-j${N}" "-R" "^benchmark" "-C" "Debug")
else ()
    add_custom_target (tests      COMMAND ${CMAKE_CTEST_COMMAND} "-j${N}" "-R" "^test")
    add_custom_target (partests   COMMAND ${CMAKE_CTEST_COMMAND} "-R" "^partest")
    add_custom_target (benchmarks COMMAND ${CMAKE_CTEST_COMMAND} "-j${N}" "-R" "^benchmark")
endif ()
endif(USE_OOFEM_EXE)

######################## Tests of Python bindings ########################################
if (USE_PYBIND_BINDINGS)
    add_test (NAME "bindings/python/tests" WORKING_DIRECTORY ${oofem_SOURCE_DIR}/bindings/python/tests COMMAND ${Python3_EXECUTABLE} "-m" "pytest" ${oofem_SOURCE_DIR}/bindings/python/tests)
    set_tests_properties("bindings/python/tests" PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}:${oofem_SOURCE_DIR}/bindings/python/oofem)
endif()
######################## End of tests of Python bindings ########################################