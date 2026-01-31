# tests/regression/CMakeLists.txt

# Ensure CTest is enabled
include(CTest)

# Find Python if some regression tests use Python scripts
find_package(Python3 COMPONENTS Interpreter)

if(USE_OOFEM_EXE)
set (oofem_REGRESSION_TEST_DIR "${oofem_SOURCE_DIR}/tests/regression")
set (oofem_BENCHMARK_DIR "${oofem_REGRESSION_TEST_DIR}/benchmark")
set (oofem_cmd $<TARGET_FILE:oofem> -Werror)
#file(TO_CMAKE_PATH "${oofem_TEST_DIR}/valgrind.supp" MEMORYCHECK_SUPPRESSIONS_FILE)
#set (CTEST_MEMORYCHECK_SUPPRESSIONS_FILE ${oofem_TEST_DIR}/tests/valgrind.supp)

if (USE_MPI_PARALLEL)
    set (par_dir ${oofem_REGRESSION_TEST_DIR}/partests)
    file (GLOB parallel_tests RELATIVE "${par_dir}" "${par_dir}/*/")
    foreach (case ${parallel_tests})
        file (GLOB files "${par_dir}/${case}/${case}.oofem.in.*")
        list (LENGTH files num_files)
        add_test (NAME "partest_${case}" WORKING_DIRECTORY ${par_dir}/${case} COMMAND "mpirun" "--oversubscribe" "-np" ${num_files} ${oofem_cmd} "-p" "-f" ${case}.oofem.in "-ksp_type" "cg")
    endforeach (case)
    # partest_brazil_2d_nl7 requires more time on testing server than default 1500
    set_tests_properties(partest_brazil_2d_nl7 PROPERTIES TIMEOUT 2500)
endif ()


function(oofem_regression_ctest dir)
    # discover tests
    file(GLOB gg RELATIVE "${oofem_REGRESSION_TEST_DIR}/${dir}" "${oofem_REGRESSION_TEST_DIR}/${dir}/*.in")
    foreach(case ${gg})
        add_test(NAME "test_${dir}/${case}" WORKING_DIRECTORY "${oofem_REGRESSION_TEST_DIR}/${dir}" COMMAND ${oofem_cmd} "-f" ${case})
    endforeach(case)

    # discover benchmarks
    file(GLOB gg RELATIVE "${oofem_BENCHMARK_DIR}/${dir}" "${oofem_BENCHMARK_DIR}/${dir}/*.in")
    foreach(case ${gg})
        add_test(NAME "benchmark_${dir}/${case}" WORKING_DIRECTORY "${oofem_BENCHMARK_DIR}/${dir}" COMMAND ${oofem_cmd} "-f" ${case})
        set_tests_properties("benchmark_${dir}/${case}" PROPERTIES LABELS "benchmark")
    endforeach(case)
    

    # repeat inside xml/ for both tests and benchmarks, if XML is enabled
    if(USE_XML)
        #file(GLOB gg RELATIVE "${oofem_TEST_DIR}/xml/${dir}" "${oofem_TEST_DIR}/xml/${dir}/*.xml")
        file(GLOB gg RELATIVE "${oofem_REGRESSION_TEST_DIR}/${dir}" "${oofem_REGRESSION_TEST_DIR}/${dir}/*.xml")
        foreach(case ${gg})
            # skip subproblems; they could be called something else than *.xml but then DataReader would not detect XML correctly
            if(${case} MATCHES "[_.]in.[a-z]+.xml")
                continue()
            endif()
            add_test(NAME "test_xml_${dir}/${case}" WORKING_DIRECTORY "${oofem_REGRESSION_TEST_DIR}/${dir}" COMMAND ${oofem_cmd} "-f" ${case})
        endforeach(case)
        file(GLOB gg RELATIVE "${oofem_BENCHMARK_DIR}/${dir}" "${oofem_BENCHMARK_DIR}/${dir}/*.xml")
        foreach(case ${gg})
            add_test(NAME "benchmark_xml_${dir}/${case}" WORKING_DIRECTORY "${oofem_BENCHMARK_DIR}/${dir}" COMMAND ${oofem_cmd} "-f" ${case})
        endforeach(case)
    endif()
endfunction()

function(add_oofem_benchmarks dir wildcard)
endfunction()


# Sequential test cases
if (USE_SM)
    oofem_regression_ctest(sm)
    oofem_regression_ctest(sm/NAFEMS)

    file (GLOB sm_tests RELATIVE "${oofem_REGRESSION_TEST_DIR}/sm" "${oofem_REGRESSION_TEST_DIR}/sm/*.sh")
    foreach (case ${sm_tests})
        add_test (NAME "test_sm_${case}" WORKING_DIRECTORY ${oofem_REGRESSION_TEST_DIR}/sm COMMAND bash ${case} ${oofem_cmd})
    endforeach (case)
endif ()


if (USE_FM)
    oofem_regression_ctest(fm)
endif()

if (USE_TM)
    oofem_regression_ctest(tm)
endif()

if (USE_AM)
    oofem_regression_ctest(am)
endif()

if (USE_TM AND USE_SM)
    oofem_regression_ctest(tmsm)
endif()

if (USE_TM AND USE_FM)
    oofem_regression_ctest(tmfm)
endif()

if (USE_TM AND USE_CEMHYD)
    oofem_regression_ctest(tmcemhyd)
endif()

if (USE_FM AND USE_SM AND USE_PFEM)
    oofem_regression_ctest(fmpfem)
endif()

if (USE_SM AND USE_MFRONT)
    oofem_regression_ctest(smmfront)
endif()

if (USE_SM AND USE_PYTHON_EXTENSION)
    oofem_regression_ctest(sm/python)
endif()

if (USE_MPM)
    oofem_regression_ctest(mpm)
endif()

if (USE_TM AND USE_SM AND USE_DSS AND USE_IML)
    oofem_regression_ctest(tmsmdssiml)
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
    add_test (NAME "test_python_pybind11_bindings" WORKING_DIRECTORY ${oofem_SOURCE_DIR}/bindings/python/tests COMMAND ${Python_EXECUTABLE} "-m" "pytest" ${oofem_SOURCE_DIR}/bindings/python/tests)
    set_tests_properties("test_python_pybind11_bindings" PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}:${oofem_SOURCE_DIR}/bindings/python/oofem)
endif()
######################## End of tests of Python bindings ########################################