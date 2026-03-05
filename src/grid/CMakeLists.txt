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

# Link with core library - this inherits Eigen includes from qlt_core
target_link_libraries(qlt_grid PRIVATE qlt_core)

# Base include directories
set(GRID_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_BINARY_DIR}/include  # For config.h
)

# Add HDF5 includes and libraries only if found
if(HDF5_FOUND)
    list(APPEND GRID_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    target_compile_definitions(qlt_grid PRIVATE -DHAVE_HDF5)
    target_link_libraries(qlt_grid PRIVATE ${HDF5_LIBRARIES})
    
    # For Windows, handle HDF5 dependencies
    if(WIN32 AND HDF5_LIBRARIES)
        target_link_libraries(qlt_grid PRIVATE 
            ${HDF5_LIBRARIES}
        )
    endif()
    
    message(STATUS "HDF5 support enabled for grid module")
else()
    message(STATUS "Building grid module without HDF5 support")
endif()

# Add _USE_MATH_DEFINES for Windows to get M_PI
target_compile_definitions(qlt_grid PRIVATE _USE_MATH_DEFINES)

# Include directories
target_include_directories(qlt_grid PRIVATE ${GRID_INCLUDE_DIRS})

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