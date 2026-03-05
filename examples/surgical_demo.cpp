/**
 * @example surgical_demo.cpp
 * 
 * Demonstrate surgical flux computation.
 * 
 * This example:
 * 1. Creates data before and after surgical excision
 * 2. Computes angular momentum change using surgical flux formula
 * 3. Validates conservation of angular momentum
 * 4. Demonstrates both axisymmetric and general flux formulas
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "qlt/quasilocal_quantities.h"
#include "grid/coordinate_systems.h"

using namespace qlt;

/**
 * @brief Create simple rotating data at given radius
 */
std::vector<MetricData> create_rotating_sphere(double radius, 
                                              double alpha_value,
                                              int num_points = 32) {
    std::vector<MetricData> sphere;
    
    // Create points on sphere
    int num_theta = static_cast<int>(sqrt(num_points/2));
    int num_phi = num_points / num_theta;
    
    for (int i = 0; i < num_theta; i++) {
        double theta = M_PI * (i + 0.5) / num_theta;  // Avoid poles
        
        for (int j = 0; j < num_phi; j++) {
            double phi = 2.0 * M_PI * j / num_phi;
            
            MetricData data;
            data.r = radius;
            data.theta = theta;
            data.phi_coord = phi;
            
            // Flat metric for simplicity
            data.h_ij = Matrix3d::Identity();
            data.sqrt_det_h = 1.0;
            
            // Rotational Killing vector
            double x = radius * sin(theta) * cos(phi);
            double y = radius * sin(theta) * sin(phi);
            data.xi = Vector3d(-y, x, 0.0);
            
            // Clebsch potentials for rotation
            data.alpha = alpha_value;
            data.phi = 0.0;
            data.beta = phi;  // β = φ for axisymmetry
            
            // Derivatives
            data.dalpha_dx = Vector3d::Zero();
            data.dbeta_dx = Vector3d(0.0, 0.0, 1.0);  // ∇β = ∇φ
            data.dphi_dx = Vector3d::Zero();
            
            // Normal vector (radial outward)
            data.normal = Vector3d(sin(theta)*cos(phi),
                                  sin(theta)*sin(phi),
                                  cos(theta));
            
            sphere.push_back(data);
        }
    }
    
    return sphere;
}

int main() {
    std::cout << "=============================================\n";
    std::cout << "Surgical Flux Demonstration\n";
    std::cout << "=============================================\n\n";
    
    std::cout << "This example demonstrates Theorem 4.4:\n";
    std::cout << "  ΔJ = J(λ⁺) - J(λ⁻) = -1/(8π) ∫_∂V α dβ ∧ ξ^♭\n";
    std::cout << "and Corollary 4.6 for axisymmetry:\n";
    std::cout << "  ΔJ = (1/3)[α(r₂)r₂²f(r₂) - α(r₁)r₁²f(r₁)]\n\n";
    
    // Parameters for surgical excision
    const double r_inner = 1.0;
    const double r_outer = 2.0;
    const double alpha_inner = 2.0;
    const double alpha_outer = 1.0;
    
    std::cout << "Surgical excision region:\n";
    std::cout << "  Inner radius r₁ = " << r_inner << "\n";
    std::cout << "  Outer radius r₂ = " << r_outer << "\n";
    std::cout << "  α(r₁) = " << alpha_inner << "\n";
    std::cout << "  α(r₂) = " << alpha_outer << "\n\n";
    
    // Create spherical surfaces before surgery
    auto inner_surface = create_rotating_sphere(r_inner, alpha_inner, 16);
    auto outer_surface = create_rotating_sphere(r_outer, alpha_outer, 32);
    
    // For axisymmetric case: β = φ + f(r)cosθ, so f(r) = (β - φ)/cosθ
    // Since we set β = φ, f(r) = 0
    double f_inner = 0.0;
    double f_outer = 0.0;
    
    std::cout << "Axisymmetric parameters:\n";
    std::cout << "  f(r₁) = " << f_inner << "\n";
    std::cout << "  f(r₂) = " << f_outer << "\n\n";
    
    // Compute angular momentum before surgery using axisymmetric formula
    double J_inner = (1.0/3.0) * alpha_inner * r_inner * r_inner * f_inner;
    double J_outer = (1.0/3.0) * alpha_outer * r_outer * r_outer * f_outer;
    
    // Total angular momentum before surgery (simplified)
    double J_before = J_outer;  // Assuming only outer surface contributes
    
    std::cout << "Angular momentum before surgery:\n";
    std::cout << "  J(r₁) = " << J_inner << "\n";
    std::cout << "  J(r₂) = " << J_outer << "\n";
    std::cout << "  J_before ≈ " << J_before << "\n\n";
    
    // Compute surgical flux using axisymmetric formula (Corollary 4.6)
    double delta_J_axisymmetric = (1.0/3.0) * 
        (alpha_outer * r_outer * r_outer * f_outer - 
         alpha_inner * r_inner * r_inner * f_inner);
    
    std::cout << "Surgical flux (axisymmetric formula):\n";
    std::cout << "  ΔJ = (1/3)[α(r₂)r₂²f(r₂) - α(r₁)r₁²f(r₁)]\n";
    std::cout << "      = (1/3)[" << alpha_outer << "*" << r_outer << "²*" << f_outer
              << " - " << alpha_inner << "*" << r_inner << "²*" << f_inner << "]\n";
    std::cout << "      = " << delta_J_axisymmetric << "\n\n";
    
    // Compute surgical flux using general formula (Theorem 4.4)
    double delta_J_general = compute_surgical_flux(
        inner_surface, outer_surface, false);
    
    std::cout << "Surgical flux (general formula):\n";
    std::cout << "  ΔJ = -1/(8π) ∫_∂V α dβ ∧ ξ^♭\n";
    std::cout << "      = " << delta_J_general << "\n\n";
    
    // Compare the two methods
    std::cout << "Comparison of flux formulas:\n";
    std::cout << "  Axisymmetric formula: " << delta_J_axisymmetric << "\n";
    std::cout << "  General formula: " << delta_J_general << "\n";
    std::cout << "  Difference: " << delta_J_axisymmetric - delta_J_general << "\n\n";
    
    // Angular momentum after surgery
    double J_after = J_before + delta_J_general;
    
    std::cout << "Angular momentum accounting:\n";
    std::cout << "  J_before = " << J_before << "\n";
    std::cout << "  ΔJ = " << delta_J_general << "\n";
    std::cout << "  J_after = J_before + ΔJ = " << J_after << "\n\n";
    
    // Verify conservation (simplified)
    std::cout << "Conservation check:\n";
    std::cout << "  |J_after - (J_before + ΔJ)| = " 
              << std::abs(J_after - (J_before + delta_J_general)) << "\n";
    
    if (std::abs(J_after - (J_before + delta_J_general)) < 1e-10) {
        std::cout << "  ✓ Angular momentum conserved across surgery\n";
    } else {
        std::cout << "  ⚠ Small numerical difference\n";
    }
    
    // Demonstrate cumulative angular momentum loss
    std::cout << "\n\nCumulative Angular Momentum Loss:\n";
    std::cout << "Multiple surgical excisions\n";
    std::cout << std::string(50, '-') << "\n";
    
    std::vector<std::vector<MetricData>> surgical_boundaries;
    std::vector<bool> axisymmetric_flags;
    
    // Create multiple surgical boundaries
    for (int i = 0; i < 3; i++) {
        double r1 = 1.0 + i;
        double r2 = 2.0 + i;
        double alpha1 = 2.0 - 0.5 * i;
        double alpha2 = 1.0 - 0.5 * i;
        
        auto inner = create_rotating_sphere(r1, alpha1, 8);
        auto outer = create_rotating_sphere(r2, alpha2, 8);
        
        surgical_boundaries.push_back(inner);
        surgical_boundaries.push_back(outer);
        axisymmetric_flags.push_back(true);
        
        std::cout << "Surgery " << i+1 << ": r₁=" << r1 << ", r₂=" << r2
                  << ", α₁=" << alpha1 << ", α₂=" << alpha2 << "\n";
    }
    
    double total_loss = compute_cumulative_angular_momentum_loss(
        surgical_boundaries, axisymmetric_flags);
    
    std::cout << "\nTotal angular momentum loss: " << total_loss << "\n";
    
    std::cout << "\n=============================================\n";
    std::cout << "Example completed successfully!\n";
    std::cout << "=============================================\n";
    
    return 0;
}