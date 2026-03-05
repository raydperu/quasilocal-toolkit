# Example programs for quasilocal toolkit
option(BUILD_EXAMPLES "Build example programs" ON)

if(BUILD_EXAMPLES)
    # Kerr validation example
    add_executable(example_kerr kerr_validation.cpp)
    target_link_libraries(example_kerr PRIVATE qlt_core qlt_grid)
    
    # Schwarzschild mass example
    add_executable(example_schwarzschild schwarzschild_mass.cpp)
    target_link_libraries(example_schwarzschild PRIVATE qlt_core qlt_grid)
    
    # Surgical flux example
    add_executable(example_surgical surgical_demo.cpp)
    target_link_libraries(example_surgical PRIVATE qlt_core qlt_grid)
    
    # Installation
    install(TARGETS
        example_kerr
        example_schwarzschild
        example_surgical
        DESTINATION bin
    )
    
    # Install Python example
    install(FILES python_demo.py
        DESTINATION share/quasilocal-toolkit/examples
    )
    
    message(STATUS "Examples will be built")
endif()