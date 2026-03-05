# Python bindings module

# Find pybind11
find_package(pybind11 REQUIRED)

# Create the Python module
pybind11_add_module(pyqlt
    qlt_bindings.cpp
)

# Link with core and grid libraries
target_link_libraries(pyqlt PRIVATE 
    qlt_core 
    qlt_grid
    pybind11::headers
)

# Include directories
target_include_directories(pyqlt PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_BINARY_DIR}/include  # For config.h
)

# Add HDF5 support if available
if(HDF5_FOUND)
    target_include_directories(pyqlt PRIVATE ${HDF5_INCLUDE_DIRS})
    target_compile_definitions(pyqlt PRIVATE -DHAVE_HDF5)
    target_link_libraries(pyqlt PRIVATE ${HDF5_LIBRARIES})
endif()

# Add _USE_MATH_DEFINES for Windows to get M_PI
target_compile_definitions(pyqlt PRIVATE _USE_MATH_DEFINES)

# Set Python module properties
set_target_properties(pyqlt PROPERTIES
    PREFIX ""
    SUFFIX ".pyd"
    OUTPUT_NAME "qlt"
    CXX_VISIBILITY_PRESET "hidden"
    VISIBILITY_INLINES_HIDDEN 1
)

# Copy to build directory for easy testing
add_custom_command(TARGET pyqlt POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_FILE:pyqlt> 
        ${CMAKE_CURRENT_BINARY_DIR}/qlt.pyd
    COMMENT "Copying Python module to build directory"
)