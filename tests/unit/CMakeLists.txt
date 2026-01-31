add_executable(test_mpmevaluator test_mpmevaluator.C)
target_link_libraries(test_mpmevaluator PRIVATE Catch2::Catch2WithMain liboofem)

include(CTest) 
include(Catch)
catch_discover_tests(test_mpmevaluator)
#set_tests_properties(test_mpmevaluator PROPERTIES LABELS "unit")