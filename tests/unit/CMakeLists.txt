add_executable(test_mpmevaluator test_mpmevaluator.C)
target_link_libraries(test_mpmevaluator PRIVATE Catch2::Catch2WithMain liboofem)

include(CTest) 
include(Catch)
catch_discover_tests(test_mpmevaluator PROPERTIES LABELS "unit" PROPERTIES ENVIRONMENT "PATH=$<TARGET_FILE_DIR:oofem>;$<TARGET_FILE_DIR:liboofem>;$ENV{PATH}")
#set_tests_properties(test_mpmevaluator.C PROPERTIES LABELS "unit")