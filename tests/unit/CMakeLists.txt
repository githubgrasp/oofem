add_executable(test_mpmevaluator test_mpmevaluator.C)
target_link_libraries(test_mpmevaluator PRIVATE Catch2::Catch2WithMain liboofem)

include(CTest) 
include(Catch)
set(CATCH_DISCOVER_TESTS_ENVIRONMENT "PATH=$<TARGET_FILE_DIR:oofem>;$ENV{PATH}" )
catch_discover_tests(test_mpmevaluator PROPERTIES LABELS "unit")


