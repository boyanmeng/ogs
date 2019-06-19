include(CodeCoverage)
append_coverage_compiler_flags()

set(COVERAGE_GCOVR_EXCLUDES
    ${PROJECT_BINARY_DIR}/.*
    Applications/CLI/.*
    ProcessLib/.*
    .*Tests/.*
    ThirdParty/.*
)

if(LCOV_PATH AND GENHTML_PATH)
    setup_target_for_coverage_lcov(
        NAME testrunner_coverage
        EXECUTABLE ${CMAKE_COMMAND} --build . --target tests
    )
    setup_target_for_coverage_lcov(
        NAME ctest_coverage
        EXECUTABLE ${CMAKE_COMMAND} --build . --target ctest-serial
    )
else()
    message(STATUS "No lcov coverage report generated because lcov or genhtml was not found.")
endif()

if(Python_EXECUTABLE)
    setup_target_for_coverage_gcovr_xml(
        NAME testrunner_coverage_cobertura
        EXECUTABLE ${CMAKE_COMMAND} --build . --target tests
    )
    setup_target_for_coverage_gcovr_xml(
        NAME ctest_coverage_cobertura
        EXECUTABLE ${CMAKE_COMMAND} --build . --target ctest-serial
    )
else()
    message(STATUS "No cobertura coverage report generated because Python executable was not found.")
endif()

if(UNIX)
    add_custom_target(clean_coverage find . -name '*.gcda' -delete)
endif()
