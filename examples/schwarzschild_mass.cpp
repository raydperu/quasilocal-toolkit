/**
 * @example schwarzschild_mass.cpp
 * 
 * Compute quasilocal mass for Schwarzschild metric.
 * 
 * This example:
 * 1. Creates Schwarzschild metric data at various radii
 * 2. Computes quasilocal mass using Brown-York formalism
 * 3. Compares with expected mass M
 * 4. Demonstrates mass aspect function
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "qlt/quasilocal_quantities.h"
#include "grid/coordinate_systems.h"

using namespace qlt;

/**
 * @brief Create Schwarzschild metric data at given radius
 */
MetricData create_schwarzschild_data_at_radius(double r, double theta, double M) {
    MetricData data;
    
    data.r = r;
    data.theta = theta;
    data.phi_coord = 0.0;
    
    // Schwarzschild metric in Schwarzschild coordinates
    // g_rr = 1/(1 - 2M/r), g_θθ = r², g_φφ = r² sin²θ
    double grr = 1.0 / (1.0 - 2.0*M/r);
    
    data.h_ij = Matrix3d::Zero();
    data.h_ij(0,0) = grr;                         // g_rr
    data.h_ij(1,1) = r * r;                       // g_θθ
    data.h_ij(2,2) = r * r * sin(theta)*sin(theta); // g_φφ
    
    data.sqrt_det_h = sqrt(grr) * r * r * sin(theta);
    
    // Normal vector (radial outward)
    data.normal = Vector3d(1.0, 0.0, 0.0);
    
    // For Schwarzschild, angular momentum is zero
    // This implies either α = 0 or dβ = 0 (or both)
    data.alpha = 0.0;  // Choose α = 0
    data.phi = 0.0;
    data.beta = data.phi_coord;  // β = φ
    
    return data;
}

int main() {
    std::cout << "=============================================\n";
    std::cout << "Schwarzschild Mass Computation Example\n";
    std::cout << "=============================================\n\n";
    
    // Parameters
    const double M = 1.0;  // Schwarzschild mass
    
    std::cout << "Parameters:\n";
    std::cout << "  Schwarzschild mass M = " << M << "\n\n";
    
    // Create spherical surfaces at various radii
    std::vector<double> radii = {3.0, 5.0, 10.0, 20.0, 50.0, 100.0};
    const int num_theta = 4;
    const int num_phi = 8;
    
    std::cout << std::setw(10) << "Radius" 
              << std::setw(15) << "r/M" 
              << std::setw(20) << "Quasilocal Mass" 
              << std::setw(15) << "Error" 
              << std::setw(15) << "Rel. Error\n";
    std::cout << std::string(75, '-') << "\n";
    
    for (double r : radii) {
        // Create spherical surface data
        std::vector<MetricData> surface_data;
        
        for (int i = 0; i < num_theta; i++) {
            double theta = M_PI * (i + 0.5) / num_theta;  // Avoid poles
            
            for (int j = 0; j < num_phi; j++) {
                double phi = 2.0 * M_PI * j / num_phi;
                
                auto data = create_schwarzschild_data_at_radius(r, theta, M);
                data.phi_coord = phi;
                
                // Update normal vector for actual position
                data.normal = Vector3d(sin(theta)*cos(phi),
                                      sin(theta)*sin(phi),
                                      cos(theta));
                
                surface_data.push_back(data);
            }
        }
        
        // Compute quasilocal mass
        double mass = compute_quasilocal_mass(surface_data, true);
        
        // Error
        double error = mass - M;
        double rel_error = std::abs(error) / M;
        
        std::cout << std::setw(10) << std::setprecision(1) << std::fixed << r
                  << std::setw(15) << r/M
                  << std::setw(20) << std::setprecision(6) << mass
                  << std::setw(15) << error
                  << std::setw(15) << rel_error << "\n";
    }
    
    // Demonstrate mass aspect function
    std::cout << "\n\nMass Aspect Function M(r):\n";
    std::cout << "For spherical symmetry: M(r) = (r/2)(1 - 1/g_rr)\n";
    std::cout << std::string(50, '-') << "\n";
    
    std::vector<MetricData> surface_data_at_radii;
    for (double r : radii) {
        auto data = create_schwarzschild_data_at_radius(r, M_PI/2.0, M);
        surface_data_at_radii.push_back(data);
    }
    
    auto mass_aspect = compute_mass_aspect(surface_data_at_radii);
    
    std::cout << std::setw(10) << "Radius" 
              << std::setw(15) << "M(r) computed" 
              << std::setw(15) << "Expected M\n";
    std::cout << std::string(40, '-') << "\n";
    
    for (size_t i = 0; i < radii.size(); i++) {
        std::cout << std::setw(10) << std::setprecision(1) << radii[i]
                  << std::setw(15) << std::setprecision(6) << mass_aspect[i]
                  << std::setw(15) << M << "\n";
    }
    
    // ADM mass approximation at large radius
    std::cout << "\nADM Mass Approximation at Large Radius:\n";
    
    double r_large = 1000.0 * M;
    auto data_large = create_schwarzschild_data_at_radius(r_large, M_PI/2.0, M);
    std::vector<MetricData> large_surface = {data_large};
    
    double adm_mass = compute_adm_mass_approximation(large_surface);
    
    std::cout << "  At r = " << r_large << " = " << r_large/M << " M\n";
    std::cout << "  ADM mass approximation = " << adm_mass << "\n";
    std::cout << "  Error = " << adm_mass - M << "\n";
    std::cout << "  Relative error = " << std::abs(adm_mass - M)/M << "\n";
    
    std::cout << "\n=============================================\n";
    std::cout << "Example completed successfully!\n";
    std::cout << "=============================================\n";
    
    return 0;
}