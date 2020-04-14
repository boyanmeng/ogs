append_source_files(SOURCES)

ogs_add_library(ThermalTwoPhaseFlowComponential ${SOURCES})
target_link_libraries(ThermalTwoPhaseFlowComponential
                      PUBLIC ProcessLib
                      PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
