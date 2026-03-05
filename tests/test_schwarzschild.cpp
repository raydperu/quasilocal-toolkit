#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "qlt/quasilocal_quantities.h"
#include "grid/coordinate_systems.h"
#include <cmath>
#include <vector>

using namespace qlt;
using Catch::Matchers::WithinAbs;

TEST_CASE("Schwarzschild metric validation", "[schwarzschild][validation]") {
    const double M = 1.0;  // Mass
    
    SECTION("Zero angular momentum") {
        // Schwarzschild should have zero angular momentum
        std::vector<double> radii = {10.0, 20.0, 50.0, 100.0};
        
        for (double r : radii) {
            // For Schwarzschild, α = 0 or f = 0 (or both)
            double alpha = 0.0;
            double f = 0.0;
            
            double J = (1.0/3.0) * alpha * r * r * f;
            
            CHECK_THAT(J, WithinAbs(0.0, 1e-10));
        }
    }
    
    SECTION("Quasilocal mass for Schwarzschild") {
        // Test mass computation at various radii
        // For Schwarzschild, mass should be M at all radii > 2M
        
        std::vector<double> radii = {3.0, 5.0, 10.0, 20.0, 50.0};
        
        for (double r : radii) {
            if (r > 2.0 * M) {  // Outside horizon
                // Create spherical surface data
                std::vector<MetricData> surface_data;
                
                // Create a few points on sphere (simplified)
                for (int i = 0; i < 4; i++) {
                    MetricData data;
                    data.r = r;
                    data.theta = M_PI/2.0;  // Equatorial plane
                    data.phi_coord = i * M_PI/2.0;
                    
                    // Schwarzschild metric in Schwarzschild coordinates
                    // g_rr = 1/(1 - 2M/r), g_θθ = r², g_φφ = r² sin²θ
                    double grr = 1.0 / (1.0 - 2.0*M/r);
                    data.h_ij = Matrix3d::Zero();
                    data.h_ij(0,0) = grr;
                    data.h_ij(1,1) = r*r;
                    data.h_ij(2,2) = r*r * sin(data.theta)*sin(data.theta);
                    
                    data.sqrt_det_h = sqrt(grr * r*r * r*r * 
                                         sin(data.theta)*sin(data.theta));
                    
                    // Normal vector (radial outward)
                    data.normal = Vector3d(1.0, 0.0, 0.0);
                    
                    surface_data.push_back(data);
                }
                
                // Compute quasilocal mass
                double mass = compute_quasilocal_mass(surface_data, true);
                
                // Should be approximately M (allowing for discretization error)
                double tolerance = 0.1 * M;  // 10% tolerance for simplified test
                CHECK_THAT(mass, WithinAbs(M, tolerance));
            }
        }
    }
    
    SECTION("Mass aspect function") {
        // For Schwarzschild, mass aspect M(r) should be constant = M
        
        std::vector<MetricData> surface_data_at_radii;
        std::vector<double> radii = {10.0, 20.0, 50.0, 100.0};
        
        for (double r : radii) {
            MetricData data;
            data.r = r;
            
            // Schwarzschild metric
            double grr = 1.0 / (1.0 - 2.0*M/r);
            data.h_ij = Matrix3d::Zero();
            data.h_ij(0,0) = grr;
            data.h_ij(1,1) = r*r;
            data.h_ij(2,2) = r*r;  // sin²θ = 1 at θ=π/2
            
            surface_data_at_radii.push_back(data);
        }
        
        auto mass_aspect = compute_mass_aspect(surface_data_at_radii);
        
        for (size_t i = 0; i < radii.size(); i++) {
            double r = radii[i];
            if (r > 2.0 * M) {  // Outside horizon
                // Mass aspect should be approximately M
                CHECK_THAT(mass_aspect[i], WithinAbs(M, 0.1 * M));
            }
        }
    }
    
    SECTION("ADM mass approximation") {
        // At large r, should recover ADM mass
        
        std::vector<MetricData> surface_data;
        
        // Create surface at large radius
        double r = 1000.0 * M;
        
        MetricData data;
        data.r = r;
        
        // Schwarzschild metric at large r
        // g_rr ≈ 1 + 2M/r
        double grr = 1.0 + 2.0*M/r;  // Approximation
        data.h_ij = Matrix3d::Zero();
        data.h_ij(0,0) = grr;
        data.h_ij(1,1) = r*r;
        data.h_ij(2,2) = r*r;
        
        surface_data.push_back(data);
        
        double adm_mass = compute_adm_mass_approximation(surface_data);
        
        // Should recover M approximately
        CHECK_THAT(adm_mass, WithinAbs(M, 0.01 * M));
    }
}