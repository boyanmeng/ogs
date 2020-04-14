append_source_files(SOURCES)

ogs_add_library(ThermalTwoPhaseFlowComponential ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS ThermalTwoPhaseFlowComponential
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(ThermalTwoPhaseFlowComponential
                      PUBLIC ProcessLib
                      PRIVATE ParameterLib)

if(BUILD_TESTING)
    include(Tests.cmake)
endif()
