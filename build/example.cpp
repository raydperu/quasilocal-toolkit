#include <iostream>
#include <vector>
#include <cmath>
#include "axisymmetric_j.h"
#include "quasilocal_mass.h"

int main() {
    std::cout << "=== QuasiLocal Toolkit Example ===\n\n";
    
    // Example 1: Axisymmetric Angular Momentum
    // J(r) = (1/3) α(r) r² f(r) where f(r) = β - φ
    std::cout << "1. Computing Axisymmetric Angular Momentum J(r)\n";
    
    // Create test data for a Kerr-like spacetime
    std::vector<double> r = {1.0, 2.0, 3.0, 4.0, 5.0};  // Radial coordinates
    std::vector<double> alpha = {0.9, 0.95, 0.98, 0.99, 1.0};  // Lapse function
    std::vector<double> f = {0.1, 0.2, 0.3, 0.4, 0.5};  // Phase function f(r) = β - φ
    
    std::cout << "Input data:\n";
    for (size_t i = 0; i < r.size(); i++) {
        std::cout << "  r = " << r[i] << ", α = " << alpha[i] << ", f = " << f[i] << "\n";
    }
    
    // Compute J(r)
    std::vector<double> J = qlt::compute_axisymmetric_j(alpha, f, r);
    
    std::cout << "\nResults J(r):\n";
    for (size_t i = 0; i < J.size(); i++) {
        std::cout << "  J(" << r[i] << ") = " << J[i] << "\n";
    }
    
    // Example 2: Surgical Flux between two surfaces
    std::cout << "\n2. Computing Surgical Flux ΔJ between r=2 and r=4\n";
    double J_flux = qlt::compute_axisymmetric_surgical_flux(
        alpha[1], r[1], f[1],  // Surface at r=2
        alpha[3], r[3], f[3]   // Surface at r=4
    );
    std::cout << "  ΔJ = " << J_flux << "\n";
    
    // Example 3: Quasilocal Mass (simplified)
    std::cout << "\n3. Computing Quasilocal Mass\n";
    // Note: This requires MetricData structures, which need proper initialization
    // This is a placeholder to show the API
    
    std::cout << "\n=== Example Complete ===\n";
    return 0;
}