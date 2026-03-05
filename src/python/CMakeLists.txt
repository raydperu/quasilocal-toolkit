# Python bindings for quasilocal toolkit
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

if(BUILD_PYTHON_BINDINGS)
    # Find Python and PyBind11
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    find_package(pybind11 REQUIRED)
    
    message(STATUS "Python3 found: ${Python3_EXECUTABLE}")
    message(STATUS "Python3 include: ${Python3_INCLUDE_DIRS}")
    message(STATUS "Python3 library: ${Python3_LIBRARIES}")
    
    # PyBind11 module
    pybind11_add_module(pyqlt
        bindings.cpp
        python_helpers.cpp
    )
    
    # Link with core libraries
    target_link_libraries(pyqlt PRIVATE
        qlt_core
        qlt_grid
    )
    
    # Include directories
    target_include_directories(pyqlt PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../include
        ${Python3_INCLUDE_DIRS}
        ${pybind11_INCLUDE_DIRS}
    )
    
    # Set output name and properties
    set_target_properties(pyqlt PROPERTIES
        OUTPUT_NAME "qlt"
        SUFFIX "${Python3_MODULE_EXTENSION}"
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN 1
    )
    
    # Copy __init__.py to build directory
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in
        ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
        COPYONLY
    )
    
    # Install Python module
    install(TARGETS pyqlt
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/qlt
    )
    
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/qlt
    )
    
    # Python package metadata
    set(PYTHON_PACKAGE_NAME "quasilocal-toolkit")
    set(PYTHON_PACKAGE_VERSION "${PROJECT_VERSION}")
    
    # Create setup.py for pip installation
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in
        ${CMAKE_CURRENT_BINARY_DIR}/setup.py
    )
    
    # Create pip package target
    add_custom_target(pip_package
        COMMAND ${Python3_EXECUTABLE} setup.py bdist_wheel
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building Python wheel package"
    )
    
    message(STATUS "Python bindings will be built as module: qlt")
else()
    message(STATUS "Python bindings disabled")
endif()