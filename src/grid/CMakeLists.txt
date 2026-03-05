# Grid module for data I/O and geometric computations
set(GRID_SOURCES
    metric_reader.cpp
    coordinate_systems.cpp
    killing_vector.cpp
    surface_integral.cpp
    volume_integral.cpp
)

set(GRID_HEADERS
    metric_reader.h
    coordinate_systems.h
    killing_vector.h
    surface_integral.h
    volume_integral.h
)

# Create grid library
add_library(qlt_grid STATIC ${GRID_SOURCES})

# Link with core library
target_link_libraries(qlt_grid PRIVATE qlt_core)

# Include directories
target_include_directories(qlt_grid PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${EIGEN3_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

target_include_directories(qlt_grid PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# Compile options
target_compile_features(qlt_grid PUBLIC cxx_std_17)

# Installation
install(TARGETS qlt_grid
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(FILES ${GRID_HEADERS}
    DESTINATION include/qlt/grid
)