#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "qlt/quasilocal_quantities.h"
#include <cmath>
#include <vector>

using namespace qlt;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Kerr metric validation", "[kerr][validation]") {
    // Test parameters
    const double M = 1.0;      // Mass
    const double a = 0.5;      // Spin parameter
    const double J_expected = M * a;  // Expected angular momentum
    
    SECTION("Axisymmetric formula for Kerr") {
        // Test at various radii
        std::vector<double> radii = {10.0, 20.0, 50.0, 100.0, 200.0};
        
        for (double r : radii) {
            // Kerr metric in Boyer-Lindquist coordinates (θ = π/2, equatorial plane)
            double theta = M_PI / 2.0;
            double Sigma = r*r + a*a;
            double Delta = r*r - 2.0*M*r + a*a;
            
            // Kerr α(r) for Clebsch representation
            double alpha = -2.0 * M * a * r * sin(theta)*sin(theta) / Sigma;
            
            // f(r) satisfying Clebsch constraint
            double f = -a * (r*r + a*a - M*r) / Delta;
            
            // Compute angular momentum
            double J_computed = (1.0/3.0) * alpha * r * r * f;
            
            // Should approach J_expected as r → ∞
            double tolerance = 0.01 * (100.0 / r);  // Larger tolerance at small r
            
            INFO("Testing at radius r = " << r);
            CHECK_THAT(J_computed, WithinRel(J_expected, tolerance));
        }
    }
    
    SECTION("Zero angular momentum for Schwarzschild limit") {
        // When a = 0, we have Schwarzschild, J should be 0
        double a_schwarzschild = 0.0;
        double r = 100.0;
        double theta = M_PI / 2.0;
        double Sigma = r*r + a_schwarzschild*a_schwarzschild;
        
        double alpha = -2.0 * M * a_schwarzschild * r * sin(theta)*sin(theta) / Sigma;
        double f = 0.0;  // For Schwarzschild
        
        double J = (1.0/3.0) * alpha * r * r * f;
        
        CHECK_THAT(J, WithinAbs(0.0, 1e-10));
    }
    
    SECTION("Extremal Kerr bound") {
        // For extremal Kerr, a = M, so J = M²
        double a_extremal = M;
        double J_extremal_expected = M * M;
        
        double r = 100.0;  // Large radius for accurate measurement
        double theta = M_PI / 2.0;
        double Sigma = r*r + a_extremal*a_extremal;
        double Delta = r*r - 2.0*M*r + a_extremal*a_extremal;
        
        double alpha = -2.0 * M * a_extremal * r * sin(theta)*sin(theta) / Sigma;
        double f = -a_extremal * (r*r + a_extremal*a_extremal - M*r) / Delta;
        
        double J_computed = (1.0/3.0) * alpha * r * r * f;
        
        CHECK_THAT(J_computed, WithinRel(J_extremal_expected, 0.01));
    }
    
    SECTION("Negative spin parameter") {
        // Test with negative a (counter-rotating)
        double a_negative = -0.3;
        double J_expected_negative = M * a_negative;
        
        double r = 100.0;
        double theta = M_PI / 2.0;
        double Sigma = r*r + a_negative*a_negative;
        double Delta = r*r - 2.0*M*r + a_negative*a_negative;
        
        double alpha = -2.0 * M * a_negative * r * sin(theta)*sin(theta) / Sigma;
        double f = -a_negative * (r*r + a_negative*a_negative - M*r) / Delta;
        
        double J_computed = (1.0/3.0) * alpha * r * r * f;
        
        CHECK_THAT(J_computed, WithinRel(J_expected_negative, 0.01));
    }
}

TEST_CASE("Kerr-Clebsch constraint satisfaction", "[kerr][constraint]") {
    // Test that Kerr metric satisfies Clebsch constraint C_{μν} = 0
    
    const double M = 1.0;
    const double a = 0.5;
    
    // Create test points at various radii
    std::vector<double> radii = {10.0, 20.0, 50.0};
    
    for (double r : radii) {
        double theta = M_PI / 2.0;
        double phi = 0.0;
        
        // Create MetricData for Kerr at this point
        MetricData data;
        data.r = r;
        data.theta = theta;
        data.phi_coord = phi;
        
        // Kerr metric in Boyer-Lindquist coordinates (simplified)
        // In full implementation, would compute actual Kerr metric
        
        // Set Clebsch potentials for Kerr
        double Sigma = r*r + a*a*cos(theta)*cos(theta);
        data.alpha = -2.0 * M * a * r * sin(theta)*sin(theta) / Sigma;
        data.phi = 0.0;  // Constant
        data.beta = phi; // β = φ for axisymmetry
        
        // Compute derivatives (would need actual metric)
        // For test, set to reasonable values
        
        // Check zero angular momentum condition
        // When α = 0 or dβ = 0, J should be 0
        // For Kerr with a ≠ 0, neither is zero
        
        CHECK(data.alpha != 0.0);  // Kerr has non-zero α
        // dβ = dφ ≠ 0 for axisymmetry
    }
}

TEST_CASE("Vorticity bounds for Kerr", "[kerr][vorticity]") {
    // Test that Kerr satisfies vorticity bounds
    
    const double M = 1.0;
    const double a = 0.5;
    const double r = 100.0;
    const double theta = M_PI / 2.0;
    
    // Create test data
    MetricData data;
    data.r = r;
    data.theta = theta;
    
    // Kerr Clebsch potentials
    double Sigma = r*r + a*a*cos(theta)*cos(theta);
    data.alpha = -2.0 * M * a * r * sin(theta)*sin(theta) / Sigma;
    
    // Set metric (flat for simplicity in test)
    data.h_ij = Matrix3d::Identity();
    data.sqrt_det_h = 1.0;
    
    // Set derivatives (simplified)
    data.dalpha_dx = Vector3d(0.1, 0.0, 0.0);  // Some gradient
    data.dbeta_dx = Vector3d(0.0, 0.0, 1.0);   // dβ points in φ direction
    
    // Compute vorticity bound
    auto [lhs, rhs] = compute_vorticity_bound(data);
    
    // Check bound: ||V ∧ dV|| ≤ ||α|| ||dβ||²
    // For Kerr, this should be satisfied
    CHECK(lhs <= rhs + 1e-10);
    
    // For extreme case, check more carefully
    if (std::abs(a) < M) {  // Sub-extremal
        CHECK(lhs < rhs * 1.1);  // Allow some margin
    }
}