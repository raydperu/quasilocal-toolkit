cmake_minimum_required(VERSION 3.14)
project(QuasiLocalToolkit VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Policy settings
if(POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif()

# Build options
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(BUILD_TESTS "Build test suite" ON)
option(BUILD_EXAMPLES "Build examples" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Find dependencies
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# HDF5 is optional but recommended
find_package(HDF5 QUIET COMPONENTS C CXX)
if(HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_VERSION}")
    add_definitions(-DHAVE_HDF5)
else()
    message(WARNING "HDF5 not found, some features will be disabled")
endif()

# Add subdirectories
add_subdirectory(src/core)
add_subdirectory(src/grid)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(src/python)
endif()

if(BUILD_TESTS)
    add_subdirectory(tests)
endif()

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Installation
install(DIRECTORY include/qlt DESTINATION include)

# Package configuration
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/QuasiLocalToolkitConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/QuasiLocalToolkitConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/QuasiLocalToolkitConfig.cmake
    INSTALL_DESTINATION lib/cmake/QuasiLocalToolkit
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/QuasiLocalToolkitConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/QuasiLocalToolkitConfigVersion.cmake
    DESTINATION lib/cmake/QuasiLocalToolkit
)

# Export targets
install(EXPORT QuasiLocalToolkitTargets
    FILE QuasiLocalToolkitTargets.cmake
    NAMESPACE qlt::
    DESTINATION lib/cmake/QuasiLocalToolkit
)

export(EXPORT QuasiLocalToolkitTargets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/QuasiLocalToolkitTargets.cmake
    NAMESPACE qlt::
)

# Create a config.h file
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/qlt/config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/qlt/config.h
)

# Add the binary directory to the include path
target_include_directories(qlt_core PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/include
)

# Summary
message(STATUS "=" * 50)
message(STATUS "QuasiLocal Toolkit Configuration Summary")
message(STATUS "=" * 50)
message(STATUS "Version: ${PROJECT_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "")
message(STATUS "Features:")
message(STATUS "  Python bindings: ${BUILD_PYTHON_BINDINGS}")
message(STATUS "  Tests: ${BUILD_TESTS}")
message(STATUS "  Examples: ${BUILD_EXAMPLES}")
message(STATUS "  HDF5 support: ${HDF5_FOUND}")
message(STATUS "")
message(STATUS "Libraries to be built:")
message(STATUS "  qlt_core (core mathematical functions)")
message(STATUS "  qlt_grid (data I/O and geometry)")
if(BUILD_PYTHON_BINDINGS)
    message(STATUS "  pyqlt (Python bindings)")
endif()
message(STATUS "=" * 50)