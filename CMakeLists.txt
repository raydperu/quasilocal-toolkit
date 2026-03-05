cmake_minimum_required(VERSION 3.14)
project(QuasiLocalToolkit VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add custom module path for FindEigen3.cmake and other custom modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

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
find_package(Eigen3 3.3 REQUIRED)

# HDF5 is optional but recommended - with proper handling
find_package(HDF5 QUIET COMPONENTS C CXX)
if(HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_VERSION}")
    message(STATUS "HDF5 C includes: ${HDF5_C_INCLUDE_DIRS}")
    message(STATUS "HDF5 CXX includes: ${HDF5_CXX_INCLUDE_DIRS}")
    message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")
    
    # Combine include directories
    set(HDF5_INCLUDE_DIRS ${HDF5_C_INCLUDE_DIRS} ${HDF5_CXX_INCLUDE_DIRS})
    list(REMOVE_DUPLICATES HDF5_INCLUDE_DIRS)
    
    add_definitions(-DHAVE_HDF5)
else()
    message(WARNING "HDF5 not found, some features will be disabled")
    # Create empty variables to avoid errors
    set(HDF5_INCLUDE_DIRS "" CACHE INTERNAL "HDF5 include dirs")
    set(HDF5_LIBRARIES "" CACHE INTERNAL "HDF5 libraries")
    set(HDF5_C_INCLUDE_DIRS "" CACHE INTERNAL "HDF5 C include dirs")
    set(HDF5_CXX_INCLUDE_DIRS "" CACHE INTERNAL "HDF5 CXX include dirs")
endif()

# Create an interface library for the config header
add_library(qlt_config INTERFACE)
target_include_directories(qlt_config INTERFACE 
    ${CMAKE_CURRENT_BINARY_DIR}/include
)

# Create a config.h file
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/qlt/config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/qlt/config.h
)

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

# Installation - simplified for now
install(DIRECTORY include/qlt DESTINATION include)

# Package configuration - COMMENTED OUT FOR NOW
# include(CMakePackageConfigHelpers)
# 
# write_basic_package_version_file(
#     ${CMAKE_CURRENT_BINARY_DIR}/QuasiLocalToolkitConfigVersion.cmake
#     VERSION ${PROJECT_VERSION}
#     COMPATIBILITY SameMajorVersion
# )
# 
# configure_package_config_file(
#     ${CMAKE_CURRENT_SOURCE_DIR}/cmake/QuasiLocalToolkitConfig.cmake.in
#     ${CMAKE_CURRENT_BINARY_DIR}/QuasiLocalToolkitConfig.cmake
#     INSTALL_DESTINATION lib/cmake/QuasiLocalToolkit
#     NO_SET_AND_CHECK_MACRO
#     NO_CHECK_REQUIRED_COMPONENTS_MACRO
# )
# 
# install(FILES
#     ${CMAKE_CURRENT_BINARY_DIR}/QuasiLocalToolkitConfig.cmake
#     ${CMAKE_CURRENT_BINARY_DIR}/QuasiLocalToolkitConfigVersion.cmake
#     DESTINATION lib/cmake/QuasiLocalToolkit
# )
# 
# # Export targets - COMMENTED OUT
# install(EXPORT QuasiLocalToolkitTargets
#     FILE QuasiLocalToolkitTargets.cmake
#     NAMESPACE qlt::
#     DESTINATION lib/cmake/QuasiLocalToolkit
# )

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
if(HDF5_FOUND)
    message(STATUS "  HDF5 version: ${HDF5_VERSION}")
endif()
message(STATUS "")
message(STATUS "Libraries to be built:")
message(STATUS "  qlt_core (core mathematical functions)")
message(STATUS "  qlt_grid (data I/O and geometry)")
if(BUILD_PYTHON_BINDINGS)
    message(STATUS "  pyqlt (Python bindings)")
endif()
message(STATUS "=" * 50)