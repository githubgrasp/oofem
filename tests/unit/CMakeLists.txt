add_executable(test_mpmevaluator test_mpmevaluator.C)
target_link_libraries(test_mpmevaluator PRIVATE Catch2::Catch2WithMain liboofem)

include(CTest) 
include(Catch)

list(JOIN LIB_DIRS ";" LIB_PATHS)

set(CATCH_DISCOVER_TESTS_ENVIRONMENT "PATH=${LIB_PATHS};$ENV{PATH}" )


if (USE_MPM)
    catch_discover_tests(test_mpmevaluator PROPERTIES LABELS "unit")
endif()


