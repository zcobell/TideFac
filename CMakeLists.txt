# ##############################################################################
# CMake Build File for TideFac
#
# Written By: Zach Cobell
#
# ##############################################################################
#
# The CMake build system enable TideFac to be deployed and built in a cross
# platform environment.
#
# ##############################################################################

# ##############################################################################
# GENERAL OPTIONS
# ##############################################################################
cmake_minimum_required(VERSION 3.0.2)
project(
        TideFac
        LANGUAGES CXX
        VERSION 0.0.1)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CheckFunctionExists)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(cmake/StandardProjectSettings.cmake)
include(cmake/PreventInSourceBuilds.cmake)

add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_14)

if (CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    option(ENABLE_BUILD_WITH_TIME_TRACE
            "Enable -ftime-trace to generate time tracing .json files on clang"
            OFF)
    if (ENABLE_BUILD_WITH_TIME_TRACE)
        add_compile_definitions(project_options INTERFACE -ftime-trace)
    endif ()
endif ()

add_library(project_warnings INTERFACE)
include(cmake/Cache.cmake)
include(cmake/CompilerWarnings.cmake)
option(TIDEFAC_ENABLE_PYTHON "Enable the Python interface" OFF)
option(TIDEFAC_ENABLE_FORTRAN "Enable the Fortran interface" OFF)
if (TIDEFAC_ENABLE_FORTRAN)
    enable_language(Fortran)
endif (TIDEFAC_ENABLE_FORTRAN)
option(TIDEFAC_PROJECT_WARNINGS "Enable compile time warnings as errors" OFF)
if (TIDEFAC_PROJECT_WARNINGS)
    set_project_warnings(project_warnings)
endif ()
#include(cmake/Sanitizers.cmake)
#enable_sanitizers(project_options)
#include(cmake/StaticAnalyzers.cmake)
# ##############################################################################

# ##############################################################################
# Enable running tests
# ##############################################################################
if (UNIX OR CYGWIN)
    enable_testing()
endif (UNIX OR CYGWIN)
# ##############################################################################

# ##############################################################################
# Enable Coverage
# ##############################################################################
# OPTION(COVERAGE "Export Code Coverage report from tests" OFF)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
if (TIDEFAC_COVERAGE)
    if (CMAKE_COMPILER_IS_GNUCXX)
        include(CodeCoverage)
        setup_target_for_coverage(tidefac_coverage tests coverage)
        set(CMAKE_CXX_FLAGS
                "${CMAKE_CXX_FLAGS} -Wall -pedantic -pthread -g -O0 -fprofile-arcs -ftest-coverage"
                )
        set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} -Wall -pedantic -pthread -g -O0 -fprofile-arcs -ftest-coverage"
                )
    endif (CMAKE_COMPILER_IS_GNUCXX)
endif (TIDEFAC_COVERAGE)
# ##############################################################################

# ##############################################################################
# Compiler flags
# ##############################################################################
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
mark_as_advanced(CLEAR CMAKE_CXX_FLAGS_RELEASE)
mark_as_advanced(CLEAR CMAKE_CXX_FLAGS_DEBUG)
mark_as_advanced(CLEAR CMAKE_C_FLAGS_RELEASE)
mark_as_advanced(CLEAR CMAKE_C_FLAGS_DEBUG)
mark_as_advanced(CLEAR CMAKE_CXX_COMPILER)
mark_as_advanced(CLEAR CMAKE_C_COMPILER)
# ##############################################################################

# ##############################################################################
# LIBRARY VERSION
# ##############################################################################
set(TIDEFAC_VERSION_MAJOR 0)
set(TIDEFAC_VERSION_MINOR 0)
set(TIDEFAC_VERSION_PATCH 1)
set(TIDEFAC_VERSION_STRING
        ${TIDEFAC_VERSION_MAJOR}.${TIDEFAC_VERSION_MINOR}.${TIDEFAC_VERSION_PATCH})
# ##############################################################################

# ##############################################################################
# CODE VERSION (GIT)
# ##############################################################################
execute_process(
        COMMAND git describe --always --tags
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_VERSION
        RESULT_VARIABLE GIT_RETURN_VALUE
        OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT "${GIT_RETURN_VALUE}" STREQUAL "0")
    set(GIT_VERSION "${TIDEFAC_VERSION_STRING}.cv")
endif ()
message(STATUS "TideFac Version: ${GIT_VERSION}")
# ##############################################################################

# ##############################################################################
# Boost
# ##############################################################################
# find_package(Boost 1.60.0) if(NOT Boost_FOUND OR PACKAGE_VERSION_UNSUITABLE)
# message(STATUS "Using included Boost Version: 1.73.0") set(Boost_INCLUDE_DIRS
# ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/boost_1_73_0) endif()
# ##############################################################################

# ##############################################################################
# TESTING
# ##############################################################################
option(TIDEFAC_BUILD_TESTS "Build test cases" OFF)
if (TIDEFAC_BUILD_TESTS)
    enable_testing()
endif ()
option(TIDEFAC_BUILD_BENCHMARK "Build benchmark" OFF)
if (TIDEFAC_BUILD_BENCHMARK)
    set(GOOGLE_BENCHMARK_HOME
            "/opt/google-benchmark"
            CACHE STRING "Location of the Google Benchmark library")
endif (TIDEFAC_BUILD_BENCHMARK)
# ##############################################################################

# ##############################################################################
# SET THE LOCATION OF TEMPORARY STATIC LIBS
# ##############################################################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles)
# ##############################################################################

# ##############################################################################
# SWIG/PYTHON/PERL Note: This will attempt to default to Python3 Give the flag
# -DPYTHON2=ON to force Python2 usage
# ##############################################################################
if (${PYTHON2})
    message(STATUS "CMake will attempt to use the Python2 interface only")
    find_package(PythonInterp 2.7 REQUIRED)
    find_package(PythonLibs 2.7 REQUIRED)
    set(PYTHONTYPE 2)
else ()
    find_package(PythonInterp 3.0)
    find_package(PythonLibs 3.0)
    if (NOT PythonLibs_FOUND AND NOT PythonInterp_FOUND)
        find_package(PythonInterp 2.7)
        find_package(PythonLibs 2.7)
        set(PYTHONTYPE 2)
    else ()
        set(PYTHONTYPE 3)
    endif ()
endif ()

find_package(SWIG 3.0)
if (SWIG_FOUND
        AND PythonLibs_FOUND
        AND PythonInterp_FOUND)
    include(${SWIG_USE_FILE})
    if (PythonLibs_FOUND)
        execute_process(
                COMMAND ${PYTHON_EXECUTABLE} -m site --user-site
                OUTPUT_VARIABLE temp_PYTHON_SITE_DIRECTORY
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(PYTHON_INSTALL_DIRECTORY
                ${temp_PYTHON_SITE_DIRECTORY}
                CACHE STRING "Location where Python libraries will be installed")
        mark_as_advanced(CLEAR PYTHON_EXECUTABLE)
    else (PythonLibs_FOUND)
        set(PYTHON_INSTALL_DIRECTORY
                "PYTHON_INSTALL_DIRECTORY-NOTFOUND"
                CACHE STRING "Location where Python libraries will be installed")
    endif (PythonLibs_FOUND)

    # FIND_PACKAGE(PerlLibs) IF(PerlLibs_FOUND) OPTION(TEST_PERL "Test the perl
    # build for SWIG. Note: Does not work yet" OFF) ENDIF(PerlLibs_FOUND)
else (
        SWIG_FOUND
        AND PythonLibs_FOUND
        AND PythonInterp_FOUND)
    message(
            WARNING "SWIG/Python libraries could not be found. PyTidefac disabled.")
endif (
        SWIG_FOUND
        AND PythonLibs_FOUND
        AND PythonInterp_FOUND)
# ##############################################################################

# ##############################################################################
# TideFac Library
# ##############################################################################
set(TIDEFAC_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/date.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/tidefac.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/tidefac_fortran.cpp)

add_library(tidefac ${TIDEFAC_SOURCES})

set(HEADER_LIST ${CMAKE_CURRENT_SOURCE_DIR}/src/tidefac.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/date.h)

target_include_directories(tidefac PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_compile_definitions(tidefac PRIVATE GIT_VERSION="${GIT_VERSION}")
target_compile_definitions(tidefac PRIVATE TIDEFAC_LIBRARY)

set_target_properties(tidefac PROPERTIES VERSION ${TIDEFAC_VERSION_STRING}
        SOVERSION ${TIDEFAC_VERSION_MAJOR})
set_target_properties(tidefac PROPERTIES PUBLIC_HEADER "${HEADER_LIST}")

if (${CMAKE_INSTALL_PREFIX} STREQUAL "/usr/local"
        OR ${CMAKE_INSTALL_PREFIX} STREQUAL "/usr/local/"
        OR ${CMAKE_INSTALL_PREFIX} STREQUAL "/usr/"
        OR ${CMAKE_INSTALL_PREFIX} STREQUAL "/usr/")
    set(HEADER_DEST "${CMAKE_INSTALL_INCLUDEDIR}/tidefac")
else ()
    set(HEADER_DEST ${CMAKE_INSTALL_INCLUDEDIR})
endif ()

write_basic_package_version_file(
        tidefacConfigVersion.cmake
        VERSION ${TIDEFAC_VERSION_STRING}
        COMPATIBILITY SameMajorVersion)
install(
        TARGETS tidefac
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT TIDEFAC_RUNTIME
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT TIDEFAC_RUNTIME
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT TIDEFAC_DEVELOPMENT
        PUBLIC_HEADER DESTINATION ${HEADER_DEST} COMPONENT TIDEFAC_DEVELOPMENT)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/tidefacConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)

set_target_properties(tidefac PROPERTIES CMAKE_CXX_VISIBILITY_PRESET hidden)
set_target_properties(tidefac PROPERTIES CMAKE_CXX_INLINES_HIDDEN YES)

target_link_libraries(tidefac project_warnings project_options)

if (APPLE)
    set(CMAKE_MACOSX_RPATH 0)
    set_target_properties(tidefac PROPERTIES INSTALL_NAME_DIR "tidefac")
    set_target_properties(tidefac PROPERTIES MACOSX_RPATH "tidefac")
endif (APPLE)

# ##############################################################################

# ##############################################################################
# TideFac Fortran
# ##############################################################################
if (TIDEFAC_ENABLE_FORTRAN)
    add_library(tidefacf ${CMAKE_CURRENT_SOURCE_DIR}/src/tidefac.F90)
    add_dependencies(tidefacf tidefac)
    target_link_libraries(tidefacf tidefac)
    set_target_properties(tidefacf PROPERTIES VERSION ${TIDEFAC_VERSION_STRING}
            SOVERSION ${TIDEFAC_VERSION_MAJOR})
    set_target_properties(tidefacf PROPERTIES PUBLIC_HEADER "${MOD_LIST}")
    set(MOD_LIST
            ${CMAKE_BINARY_DIR}/CMakeFiles/tidefacf/tidefacmodule.mod)
    set_target_properties(
            tidefacf PROPERTIES Fortran_MODULE_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/tidefacf)
    write_basic_package_version_file(
            tidefacfConfigVersion.cmake
            VERSION ${TIDEFAC_VERSION_STRING}
            COMPATIBILITY SameMajorVersion)
    install(
            TARGETS tidefacf
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT TIDEFACF_RUNTIME
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT TIDEFACF_RUNTIME
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT TIDEFACF_DEVELOPMENT
            PUBLIC_HEADER DESTINATION ${HEADER_DEST} COMPONENT TIDEFACF_DEVELOPMENT)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/tidefacfConfigVersion.cmake
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)

endif (TIDEFAC_ENABLE_FORTRAN)
# ##############################################################################

# ##############################################################################
# PyTideFac Library
# ##############################################################################
if (TIDEFAC_ENABLE_PYTHON
        AND SWIG_FOUND
        AND PythonLibs_FOUND
        AND PythonInterp_FOUND)

    if (PYTHONTYPE EQUAL 2)
        set(PYTHONFLAG "")
    else ()
        set(PYTHONFLAG "-py3")
    endif ()

    if (WIN32)

        # ...Builds a dummy cmake project Which will add some windows specific
        # include files not required by Unix/Cygwin [mingw64] Further testing
        # required to see if this is required by MSVC compilers
        file(
                WRITE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/python_wrapper.cmake
                "
            FILE(READ \${SRC} DATA)
            FILE(WRITE \${DST} \"#include <math.h>\")
            FILE(APPEND \${DST} \"\\n\\n\")
            FILE(APPEND \${DST} \"\${DATA}\")
        ")

        add_custom_command(
                OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/python_tidefac_wrap.cxx
                COMMAND
                ${SWIG_EXECUTABLE} -outdir ${CMAKE_CURRENT_BINARY_DIR} -c++ -python
                ${PYTHONFLAG} -I${CMAKE_CURRENT_SOURCE_DIR}/src -I${PYTHON_INCLUDE_PATH}
                -o ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/python_tidefac_wrap.cxx.1
                ${CMAKE_CURRENT_SOURCE_DIR}/python/tidefac.i
                COMMAND
                ${CMAKE_COMMAND}
                -DSRC=\"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/python_tidefac_wrap.cxx.1\"
                -DDST=\"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/python_tidefac_wrap.cxx\"
                -P ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/python_wrapper.cmake
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating SWIG interface to Python...")

    else (WIN32)
        add_custom_command(
                OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/python_tidefac_wrap.cxx
                COMMAND
                ${SWIG_EXECUTABLE} -outdir ${CMAKE_CURRENT_BINARY_DIR} -c++ -python
                ${PYTHONFLAG} -I${CMAKE_CURRENT_SOURCE_DIR}/src -I${PYTHON_INCLUDE_PATH}
                -o ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/python_tidefac_wrap.cxx
                ${CMAKE_CURRENT_SOURCE_DIR}/python/tidefac.i
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                COMMENT "Generating SWIG interface to Python...")
    endif (WIN32)

    add_library(PyTidefac SHARED
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/python_tidefac_wrap.cxx)
    target_include_directories(PyTidefac PRIVATE ${PYTHON_INCLUDE_PATH}
            ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(PyTidefac tidefac ${PYTHON_LIBRARY})
    set_target_properties(PyTidefac PROPERTIES PREFIX "_")
    set_target_properties(PyTidefac PROPERTIES INSTALL_NAME_DIR "PyTidefac")
    set_property(
            DIRECTORY
            APPEND
            PROPERTY ADDITIONAL_MAKE_CLEAN_FILES PyTidefac.py
            CMakeFiles/python_tidefac_wrap.cxx)

    add_dependencies(PyTidefac tidefac)

    if (WIN32)
        set_target_properties(
                PyTidefac PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        install(
                TARGETS PyTidefac
                LIBRARY DESTINATION bin
                RUNTIME DESTINATION ${PYTHON_INSTALL_DIRECTORY})
    elseif (CYGWIN)
        set_target_properties(PyTidefac PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                ${CMAKE_CURRENT_BINARY_DIR})
        install(
                TARGETS PyTidefac
                LIBRARY DESTINATION lib
                RUNTIME DESTINATION ${PYTHON_INSTALL_DIRECTORY})
    else (WIN32)
        set_target_properties(PyTidefac PROPERTIES LIBRARY_OUTPUT_DIRECTORY
                ${CMAKE_CURRENT_BINARY_DIR})
        install(TARGETS PyTidefac LIBRARY DESTINATION ${PYTHON_INSTALL_DIRECTORY})
    endif (WIN32)

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/PyTidefac.py
            DESTINATION ${PYTHON_INSTALL_DIRECTORY})

    if (WIN32)
        set_target_properties(PyTidefac PROPERTIES SUFFIX ".pyd")
        if (CMAKE_SIZEOF_VOID_P EQUAL 8)
            # 64 bits
            target_compile_definitions(PyTidefac PRIVATE "MS_WIN64")
        elseif (CMAKE_SIZEOF_VOID_P EQUAL 4)
            # 32 bits
            target_compile_definitions(PyTidefac PRIVATE "MS_WIN32")
        endif ()
    elseif (APPLE)
        set_target_properties(PyTidefac PROPERTIES SUFFIX ".so")
    endif (WIN32)
endif ()
# ##############################################################################

# ##############################################################################
# Test Cases
# ##############################################################################
if (UNIX OR CYGWIN)
    if (TIDEFAC_BUILD_TESTS)
        # ...C++ Testing
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testing)

        set(TEST_LIST "testcase.cpp" "testcase_fortran.F90")

        foreach (TESTFILE ${TEST_LIST})
            get_filename_component(TESTNAME ${TESTFILE} NAME_WE)
            add_executable(${TESTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/testing/${TESTFILE})
            add_dependencies(${TESTNAME} tidefac)
            target_include_directories(
                    ${TESTNAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
                    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/catch2)
            target_link_libraries(${TESTNAME} tidefac tidefacf)
            set_target_properties(
                    ${TESTNAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                    ${CMAKE_CURRENT_BINARY_DIR}/cxx)

            add_test(
                    NAME TEST_${TESTNAME}
                    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/cxx/${TESTNAME}
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)

            get_filename_component(EXTENSION ${TESTFILE} EXT)
            if (${EXTENSION} STREQUAL ".F90")
                set_target_properties(${TESTNAME} PROPERTIES LINKER_LANGUAGE Fortran)
                target_include_directories(${TESTNAME} PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/tidefacf)
            endif ()

            if (CYGWIN)
                set_tests_properties(
                        TEST_${TESTNAME}
                        PROPERTIES ENVIRONMENT "PATH=$ENV{PATH}:${CMAKE_CURRENT_BINARY_DIR}")
            else (CYGWIN)
                set_tests_properties(
                        TEST_${TESTNAME}
                        PROPERTIES
                        ENVIRONMENT
                        "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${CMAKE_CURRENT_BINARY_DIR}")
            endif (CYGWIN)
        endforeach ()

        # ...Python Test
        if (SWIG_FOUND)
            option(TIDEFAC_PYTHON_TESTING "Run python3 tests" OFF)
            if (TIDEFAC_PYTHON_TESTING)
                set(PYTHON_TESTNAME TEST_python3_PyTidefac)
                add_test(
                        NAME ${PYTHON_TESTNAME}
                        COMMAND ${PYTHON_EXECUTABLE}
                        ${CMAKE_SOURCE_DIR}/testing/python/adcirc_read.py
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
                set_property(
                        TEST ${PYTHON_TESTNAME}
                        PROPERTY
                        ENVIRONMENT
                        "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${CMAKE_CURRENT_BINARY_DIR}")
                set_property(
                        TEST ${PYTHON_TESTNAME}
                        APPEND
                        PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")
            endif ()
        endif (SWIG_FOUND)
    endif (TIDEFAC_BUILD_TESTS)
endif (UNIX OR CYGWIN)
# ##############################################################################

# ##############################################################################
# Benchmarking
# ##############################################################################
if (TIDEFAC_BUILD_BENCHMARK)
    add_executable(tidefac_benchmark ${CMAKE_CURRENT_SOURCE_DIR}/bench/bench.cpp)
    target_link_libraries(tidefac_benchmark tidefac benchmark pthread)
    target_link_directories(tidefac_benchmark PRIVATE
            ${GOOGLE_BENCHMARK_HOME}/lib)
    target_include_directories(
            tidefac_benchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${GOOGLE_BENCHMARK_HOME}/include)
endif (TIDEFAC_BUILD_BENCHMARK)

# ##############################################################################
