#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "qlt/quasilocal_quantities.h"
#include <cmath>

using namespace qlt;
using Catch::Matchers::WithinAbs;

TEST_CASE("Flat space tests", "[flat][validation]") {
    SECTION("Zero mass in flat space") {
        // In flat space, quasilocal mass should be zero
        
        std::vector<MetricData> surface_data;
        
        // Create spherical surface in flat space
        double r = 10.0;
        
        MetricData data;
        data.r = r;
        data.theta = M_PI/2.0;
        data.phi_coord = 0.0;
        
        // Flat metric in spherical coordinates
        data.h_ij = Matrix3d::Zero();
        data.h_ij(0,0) = 1.0;          // g_rr = 1
        data.h_ij(1,1) = r*r;          // g_θθ = r²
        data.h_ij(2,2) = r*r;          // g_φφ = r² (sin²θ = 1 at θ=π/2)
        
        data.sqrt_det_h = r*r;  // √det(g) = r² sinθ = r² at θ=π/2
        
        // Normal vector (radial outward)
        data.normal = Vector3d(1.0, 0.0, 0.0);
        
        surface_data.push_back(data);
        
        // Compute quasilocal mass with reference term (flat space subtraction)
        double mass = compute_quasilocal_mass(surface_data, true);
        
        // Should be approximately zero
        CHECK_THAT(mass, WithinAbs(0.0, 1e-10));
    }
    
    SECTION("Angular momentum for rotational Killing vector") {
        // In flat space with rotational Killing vector,
        // angular momentum should be computable
        
        std::vector<MetricData> volume_data;
        
        // Create simple volume data (single point for test)
        MetricData data;
        data.r = 5.0;
        data.theta = M_PI/2.0;
        data.phi_coord = 0.0;
        
        // Flat metric
        data.h_ij = Matrix3d::Identity();
        data.sqrt_det_h = 1.0;
        
        // Rotational Killing vector about z-axis
        double x = data.r * sin(data.theta) * cos(data.phi_coord);
        double y = data.r * sin(data.theta) * sin(data.phi_coord);
        data.xi = Vector3d(-y, x, 0.0);
        
        // Derivatives of Killing vector
        data.dxi_dx = Matrix3d::Zero();
        data.dxi_dx(0,1) = -1.0;  // ∂ξ^x/∂y = -1
        data.dxi_dx(1,0) = 1.0;   // ∂ξ^y/∂x = 1
        
        // Clebsch potentials (for flat space rotation)
        data.alpha = 1.0;  // Some amplitude
        data.phi = 0.0;
        data.beta = data.phi_coord;  // β = φ
        
        // Derivatives
        data.dalpha_dx = Vector3d::Zero();
        data.dbeta_dx = Vector3d(0.0, 0.0, 1.0);  // ∇β points in φ direction
        data.dphi_dx = Vector3d::Zero();
        
        volume_data.push_back(data);
        
        // Create empty surface data (boundary at infinity)
        std::vector<MetricData> surface_data;
        
        // Compute Clebsch-Komar angular momentum
        double J = compute_clebsch_komar_j(surface_data, volume_data, 2.0);
        
        // In flat space with our simple data, should get some value
        // Exact value depends on the volume integral
        CHECK(std::isfinite(J));
    }
    
    SECTION("Zero angular momentum condition in flat space") {
        // Test Corollary 3.3/4.3: If α ≡ 0 or dβ ≡ 0, then J = 0
        
        std::vector<MetricData> test_data;
        
        // Create test point
        MetricData data;
        data.r = 1.0;
        data.h_ij = Matrix3d::Identity();
        
        SECTION("Case 1: α = 0") {
            data.alpha = 0.0;
            data.dbeta_dx = Vector3d(1.0, 0.0, 0.0);  // Non-zero dβ
            
            test_data.push_back(data);
            
            auto [alpha_zero, dbeta_zero] = check_zero_angular_momentum(test_data, 1e-10);
            
            CHECK(alpha_zero == true);
            CHECK(dbeta_zero == false);
            
            // J should be zero
            // Would need actual J computation to verify
        }
        
        SECTION("Case 2: dβ = 0") {
            data.alpha = 1.0;  // Non-zero α
            data.dbeta_dx = Vector3d::Zero();  // dβ = 0
            
            test_data.push_back(data);
            
            auto [alpha_zero, dbeta_zero] = check_zero_angular_momentum(test_data, 1e-10);
            
            CHECK(alpha_zero == false);
            CHECK(dbeta_zero == true);
        }
        
        SECTION("Case 3: Both non-zero") {
            data.alpha = 1.0;
            data.dbeta_dx = Vector3d(1.0, 0.0, 0.0);
            
            test_data.push_back(data);
            
            auto [alpha_zero, dbeta_zero] = check_zero_angular_momentum(test_data, 1e-10);
            
            CHECK(alpha_zero == false);
            CHECK(dbeta_zero == false);
        }
    }
    
    SECTION("Vorticity bound in flat space") {
        // Test Theorem 3.1/4.1: ||V ∧ dV|| ≤ ||α|| ||dβ||²
        
        MetricData data;
        data.r = 1.0;
        data.h_ij = Matrix3d::Identity();
        data.sqrt_det_h = 1.0;
        
        // Test various values
        data.alpha = 2.0;
        data.dalpha_dx = Vector3d(0.1, 0.2, 0.3);
        data.dbeta_dx = Vector3d(1.0, 0.0, 0.0);
        data.dphi_dx = Vector3d(0.5, 0.5, 0.5);
        
        auto [lhs, rhs] = compute_vorticity_bound(data);
        
        // Check the bound
        CHECK(lhs <= rhs + 1e-10);
        
        // Also check that bound is reasonable
        CHECK(rhs >= 0.0);
        
        // Test edge case: α = 0
        MetricData data_zero_alpha = data;
        data_zero_alpha.alpha = 0.0;
        
        auto [lhs2, rhs2] = compute_vorticity_bound(data_zero_alpha);
        CHECK(rhs2 == 0.0);  ||α|| = 0, so bound is 0
        CHECK(lhs2 >= 0.0);
    }
}