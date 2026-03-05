# Test suite for quasilocal toolkit
option(BUILD_TESTS "Build test suite" ON)

if(BUILD_TESTS)
    # Find testing frameworks
    find_package(Catch2 3.0 QUIET)
    
    if(Catch2_FOUND)
        message(STATUS "Found Catch2: ${Catch2_DIR}")
        
        # Test executables
        add_executable(test_kerr test_kerr.cpp)
        add_executable(test_schwarzschild test_schwarzschild.cpp)
        add_executable(test_flat_space test_flat_space.cpp)
        add_executable(test_axisymmetric test_axisymmetric.cpp)
        add_executable(test_surgical_flux test_surgical_flux.cpp)
        add_executable(test_all test_all.cpp)
        
        # Link with libraries
        foreach(test_name test_kerr test_schwarzschild test_flat_space 
                         test_axisymmetric test_surgical_flux test_all)
            target_link_libraries(${test_name} PRIVATE
                qlt_core
                qlt_grid
                Catch2::Catch2WithMain
            )
            
            target_include_directories(${test_name} PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR}/../include
                ${EIGEN3_INCLUDE_DIRS}
            )
        endforeach()
        
        # Add tests to CTest
        add_test(NAME TestKerr COMMAND test_kerr)
        add_test(NAME TestSchwarzschild COMMAND test_schwarzschild)
        add_test(NAME TestFlatSpace COMMAND test_flat_space)
        add_test(NAME TestAxisymmetric COMMAND test_axisymmetric)
        add_test(NAME TestSurgicalFlux COMMAND test_surgical_flux)
        add_test(NAME TestAll COMMAND test_all)
        
        # Enable test discovery for Catch2
        catch_discover_tests(test_all)
        
    else()
        message(WARNING "Catch2 not found, tests will not be built")
    endif()
    
    # Add Python tests if Python bindings are built
    if(BUILD_PYTHON_BINDINGS AND Python3_EXECUTABLE)
        add_test(NAME PythonImportTest
            COMMAND ${Python3_EXECUTABLE} -c "import qlt; print('Python import OK')")
        
        add_test(NAME PythonBasicTest
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python_test.py)
    endif()
    
else()
    message(STATUS "Tests disabled")
endif()