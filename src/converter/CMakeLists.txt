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
outputmanager.C
overlap.C
prism.C
region.C
surface.C
tetra.C
vertex.C
converterinputrecord.C
convertertxtinputrecord.C
convertertxtdatareader.C
)

add_library(converter STATIC ${CONVERTER_SOURCES})
target_link_libraries(converter PUBLIC liboofem)
target_include_directories(converter PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(converter.exec main.C)
target_link_libraries(converter.exec PRIVATE converter liboofem)
