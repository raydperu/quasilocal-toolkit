#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_session.hpp>

// Include all test files
#include "test_kerr.cpp"
#include "test_schwarzschild.cpp"
#include "test_flat_space.cpp"
#include "test_axisymmetric.cpp"
#include "test_surgical_flux.cpp"

int main(int argc, char* argv[]) {
    // Global setup if needed
    
    // Run tests
    int result = Catch::Session().run(argc, argv);
    
    // Global cleanup if needed
    
    return result;
}