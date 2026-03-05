# FindEigen3.cmake - Find Eigen3 library
#
# This module defines:
#  EIGEN3_INCLUDE_DIR - where to find Eigen/Core
#  EIGEN3_FOUND - True if Eigen found
#  Eigen3::Eigen - imported target

# First check if Eigen3_DIR is set
if(Eigen3_DIR)
    set(EIGEN3_SEARCH_PATHS ${Eigen3_DIR})
else()
    set(EIGEN3_SEARCH_PATHS
        "C:/Users/PC/Downloads/eigen-3.4.0/eigen-3.4.0"
        "C:/Users/PC/Downloads/eigen-3.4.0"
        ${CMAKE_SOURCE_DIR}/eigen
        ${CMAKE_SOURCE_DIR}/third_party/eigen
        ${CMAKE_SOURCE_DIR}/external/eigen
    )
endif()

message(STATUS "Searching for Eigen3 in: ${EIGEN3_SEARCH_PATHS}")

find_path(EIGEN3_INCLUDE_DIR
    NAMES Eigen/Core
    PATHS ${EIGEN3_SEARCH_PATHS}
    NO_DEFAULT_PATH
    DOC "Eigen3 include directory"
)

if(NOT EIGEN3_INCLUDE_DIR)
    # Try with default paths as fallback
    find_path(EIGEN3_INCLUDE_DIR
        NAMES Eigen/Core
        PATH_SUFFIXES eigen3 Eigen include/eigen3 include
        DOC "Eigen3 include directory"
    )
endif()

# Check version
if(EIGEN3_INCLUDE_DIR AND EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h")
    file(STRINGS "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h" _eigen_version_lines
        REGEX "#define[ \t]+EIGEN_(WORLD|MAJOR|MINOR)_VERSION")
    
    foreach(line ${_eigen_version_lines})
        if(line MATCHES "#define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)")
            set(EIGEN3_VERSION_MAJOR "${CMAKE_MATCH_1}")
        elseif(line MATCHES "#define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)")
            set(EIGEN3_VERSION_MINOR "${CMAKE_MATCH_1}")
        elseif(line MATCHES "#define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)")
            set(EIGEN3_VERSION_PATCH "${CMAKE_MATCH_1}")
        endif()
    endforeach()
    
    if(EIGEN3_VERSION_MAJOR AND EIGEN3_VERSION_MINOR AND EIGEN3_VERSION_PATCH)
        set(EIGEN3_VERSION "${EIGEN3_VERSION_MAJOR}.${EIGEN3_VERSION_MINOR}.${EIGEN3_VERSION_PATCH}")
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen3
    REQUIRED_VARS EIGEN3_INCLUDE_DIR
    VERSION_VAR EIGEN3_VERSION
)

if(EIGEN3_FOUND AND NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}"
    )
    
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR} (version ${EIGEN3_VERSION})")
endif()

mark_as_advanced(EIGEN3_INCLUDE_DIR)