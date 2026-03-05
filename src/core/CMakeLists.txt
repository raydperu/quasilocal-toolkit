# Core quasilocal toolkit library
set(CORE_SOURCES
    axisymmetric_j.cpp
    clebsch_komar.cpp
    quasilocal_mass.cpp
    surgical_flux.cpp
    vorticity_bounds.cpp
)

set(CORE_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/qlt/quasilocal_quantities.h
    axisymmetric_j.h
    clebsch_komar.h
    quasilocal_mass.h
    surgical_flux.h
    vorticity_bounds.h
)

# Create the core library
add_library(qlt_core STATIC ${CORE_SOURCES})

# Include directories
target_include_directories(qlt_core PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${EIGEN3_INCLUDE_DIRS}
)

# Public interface
target_include_directories(qlt_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
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