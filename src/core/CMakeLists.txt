# Core quasilocal toolkit library
set(CORE_SOURCES
    axisymmetric_j.cpp
    clebsch_komar.cpp
    quasilocal_mass.cpp
    surgical_flux.cpp
    vorticity_bounds.cpp
)

set(CORE_HEADERS
    axisymmetric_j.h
    clebsch_komar.h
    quasilocal_mass.h
    surgical_flux.h
    vorticity_bounds.h
    ${CMAKE_CURRENT_SOURCE_DIR}/quasilocal_quantities.h
)

# Find Eigen3
find_package(Eigen3 REQUIRED)

# Create the core library
add_library(qlt_core STATIC ${CORE_SOURCES})

# Include directories
target_include_directories(qlt_core PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_BINARY_DIR}/include  # For config.h
)

# Link with Eigen - this automatically handles include paths
target_link_libraries(qlt_core PUBLIC Eigen3::Eigen)

# Add _USE_MATH_DEFINES for Windows to get M_PI
target_compile_definitions(qlt_core PRIVATE _USE_MATH_DEFINES)

# Public interface - export local include directory
target_include_directories(qlt_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Compile options
target_compile_features(qlt_core PUBLIC cxx_std_17)
set_target_properties(qlt_core PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN 1
)

# Installation
install(TARGETS qlt_core
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(FILES ${CORE_HEADERS}
    DESTINATION include/qlt
)