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
    generatorinputrecord.C
    generatortxtdatareader.C
    generatortxtinputrecord.C
 )

add_library(generator STATIC ${GENERATOR_SOURCES})
target_link_libraries(generator PUBLIC liboofem)
target_include_directories(generator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(generator.exec main.C)
target_link_libraries(generator.exec PRIVATE generator liboofem)

# Optional: don't build generator_exec by default with `make`
#set_target_properties(generator.exec PROPERTIES EXCLUDE_FROM_ALL TRUE)
