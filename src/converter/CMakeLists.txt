set(CONVERTER_SOURCES
curve.C
cylinder.C
ellipsoid.C
fibre.C
grid.C
gridcomponent.C
inclusion.C
interfacecylinder.C
interfacesphere.C
line.C
lognormalGenerator.C
main.C
octreegridlocalizer.C
prism.C
region.C
surface.C
vertex.C
boundarysphere.C
)

add_library(converter STATIC ${CONVERTER_SOURCES})
target_link_libraries(converter PUBLIC liboofem)
target_include_directories(converter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(converter.exec main.C)
target_link_libraries(converter.exec PRIVATE converter liboofem)

enable_testing()

# Converter tests — diff-based, one subdirectory per test case.
# Each subdirectory must contain control.in (with a `#@mesher t3d|qhull`
# directive) and reference.in. Mesh files are picked from disk:
#   mesh.t3d                                       → t3d
#   mesh.nodes + mesh.voronoi                      → qhull
# Output is always named oofem.in (standardised in test control files).

set(CONVERTER_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set(CONVERTER_RUN_TEST ${CONVERTER_TEST_DIR}/run_test.cmake)

file(GLOB _conv_test_entries RELATIVE ${CONVERTER_TEST_DIR} ${CONVERTER_TEST_DIR}/*)

foreach(_entry ${_conv_test_entries})
    set(_dir ${CONVERTER_TEST_DIR}/${_entry})
    if(NOT IS_DIRECTORY ${_dir})
        continue()
    endif()
    if(NOT EXISTS ${_dir}/control.in OR NOT EXISTS ${_dir}/reference.in)
        continue()
    endif()

    if(EXISTS ${_dir}/mesh.t3d)
        set(_mesh_args "${_dir}/mesh.t3d")
    elseif(EXISTS ${_dir}/mesh.nodes AND EXISTS ${_dir}/mesh.voronoi)
        set(_mesh_args "${_dir}/mesh.nodes;${_dir}/mesh.voronoi")
    else()
        message(WARNING "converter test '${_entry}': no recognised mesh file found, skipping")
        continue()
    endif()

    add_test(
        NAME test_converter_${_entry}
        COMMAND ${CMAKE_COMMAND}
            -DCONVERTER=$<TARGET_FILE:converter.exec>
            -DTEST_DIR=${_dir}
            "-DMESH_ARGS=${_mesh_args}"
            -DOUTPUT=oofem.in
            -P ${CONVERTER_RUN_TEST}
    )
endforeach()
