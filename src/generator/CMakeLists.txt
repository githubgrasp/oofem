set(GENERATOR_SOURCES
    grid.C
    vertex.C
    curve.C
    gridcomponent.C
    surface.C
    region.C
    gridlocalizer.C
    octreegridlocalizer.C
    inclusion.C
    interfacesphere.C
    sphere.C
    interfacecylinder.C
    cylinder.C
    interfacesurface.C
    interfaceplane.C
    prism.C
    refinement.C
    refineprism.C
    notch.C
    rect.C
   )

add_library(generator STATIC ${GENERATOR_SOURCES})
target_link_libraries(generator PUBLIC liboofem)
target_include_directories(generator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(generator.exec main.C)
target_link_libraries(generator.exec PRIVATE generator liboofem)

# Optional: don't build generator_exec by default with `make`
#set_target_properties(generator.exec PROPERTIES EXCLUDE_FROM_ALL TRUE)

enable_testing()

# Generator tests — diff-based, one subdirectory per test case.
# Each subdirectory must contain mesh.in and reference.nodes. The driver
# runs generator.exec on mesh.in and diffs the resulting nodes.dat against
# reference.nodes.

set(GENERATOR_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set(GENERATOR_RUN_TEST ${GENERATOR_TEST_DIR}/run_test.cmake)

file(GLOB _gen_test_entries RELATIVE ${GENERATOR_TEST_DIR} ${GENERATOR_TEST_DIR}/*)

foreach(_entry ${_gen_test_entries})
    set(_dir ${GENERATOR_TEST_DIR}/${_entry})
    if(NOT IS_DIRECTORY ${_dir})
        continue()
    endif()
    if(NOT EXISTS ${_dir}/mesh.in OR NOT EXISTS ${_dir}/reference.nodes)
        continue()
    endif()

    add_test(
        NAME test_generator_${_entry}
        COMMAND ${CMAKE_COMMAND}
            -DGENERATOR=$<TARGET_FILE:generator.exec>
            -DTEST_DIR=${_dir}
            -P ${GENERATOR_RUN_TEST}
    )
endforeach()
