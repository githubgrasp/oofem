include(CTest) 
include(Catch)

#list(JOIN LIB_DIRS ";" LIB_PATHS)
#message("Library dirs : ${LIB_DIRS}")
#message("Library paths for test_mpmevaluator: ${LIB_PATHS}")


#set(CATCH_DISCOVER_TESTS_ENVIRONMENT "PATH=${LIB_PATHS};$ENV{PATH}" )
if (USE_MPM)
    message (STATUS "Building unit test for MPM evaluator")
    add_executable(test_mpmevaluator test_mpmevaluator.C)
    target_link_libraries(test_mpmevaluator PRIVATE Catch2::Catch2WithMain liboofem)
    # 2. Add the COPY command to the TEST target
    add_custom_command(TARGET test_mpmevaluator POST_BUILD
        if (WIN32)
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:liboofem>
            $<TARGET_FILE_DIR:test_mpmevaluator>
            COMMENT "Copying oofem.dll to test directory to fix WIN Exit code 0xc0000135"
        endif()
    )
    catch_discover_tests(test_mpmevaluator
        PROPERTIES
            ENVIRONMENT "PATH=$<TARGET_FILE_DIR:liboofem>;$ENV{PATH}"
    )
    
    add_custom_command(TARGET test_mpmevaluator POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "DEBUG: liboofem directory is: $<TARGET_FILE_DIR:liboofem>"
)
    
endif()


