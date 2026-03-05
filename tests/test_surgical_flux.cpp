#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "qlt/quasilocal_quantities.h"
#include <cmath>
#include <vector>

using namespace qlt;
using Catch::Matchers::WithinAbs;

TEST_CASE("Surgical flux calculations", "[surgical][flux]") {
    SECTION("Axisymmetric surgical flux formula") {
        // Test Corollary 4.6: ΔJ = (1/3)[α(r₂)r₂²f(r₂) - α(r₁)r₁²f(r₁)]
        
        double alpha1 = 2.0, r1 = 1.0, f1 = 3.0;
        double alpha2 = 3.0, r2 = 2.0, f2 = 4.0;
        
        double term1 = alpha1 * r1 * r1 * f1;
        double term2 = alpha2 * r2 * r2 * f2;
        double delta_J_expected = (1.0/3.0) * (term2 - term1);
        
        // Create test surface data
        std::vector<MetricData> inner_data(1);
        std::vector<MetricData> outer_data(1);
        
        inner_data[0].alpha = alpha1;
        inner_data[0].r = r1;
        inner_data[0].beta = f1;  // Using β to store f for test
        inner_data[0].phi_coord = 0.0;
        inner_data[0].theta = M_PI/2.0;
        
        outer_data[0].alpha = alpha2;
        outer_data[0].r = r2;
        outer_data[0].beta = f2;
        outer_data[0].phi_coord = 0.0;
        outer_data[0].theta = M_PI/2.0;
        
        double delta_J = compute_surgical_flux(inner_data, outer_data, true);
        
        CHECK_THAT(delta_J, WithinAbs(delta_J_expected, 1e-10));
    }
    
    SECTION("Zero surgical flux cases") {
        SECTION("Same surface") {
            std::vector<MetricData> data(1);
            data[0].alpha = 1.0;
            data[0].r = 1.0;
            data[0].beta = 1.0;
            data[0].theta = M_PI/2.0;
            data[0].phi_coord = 0.0;
            
            double delta_J = compute_surgical_flux(data, data, true);
            
            CHECK_THAT(delta_J, WithinAbs(0.0, 1e-10));
        }
        
        SECTION("α = 0 on both surfaces") {
            std::vector<MetricData> inner_data(1);
            std::vector<MetricData> outer_data(1);
            
            inner_data[0].alpha = 0.0;
            inner_data[0].r = 1.0;
            inner_data[0].beta = 1.0;
            inner_data[0].theta = M_PI/2.0;
            
            outer_data[0].alpha = 0.0;
            outer_data[0].r = 2.0;
            outer_data[0].beta = 2.0;
            outer_data[0].theta = M_PI/2.0;
            
            double delta_J = compute_surgical_flux(inner_data, outer_data, true);
            
            CHECK_THAT(delta_J, WithinAbs(0.0, 1e-10));
        }
    }
    
    SECTION("Surgical flux sign convention") {
        // Flux should be positive when angular momentum increases outward
        
        std::vector<MetricData> inner_data(1);
        std::vector<MetricData> outer_data(1);
        
        // Case 1: α and f positive, r2 > r1
        inner_data[0].alpha = 1.0;
        inner_data[0].r = 1.0;
        inner_data[0].beta = 1.0;  // f = 1
        inner_data[0].theta = M_PI/2.0;
        inner_data[0].phi_coord = 0.0;
        
        outer_data[0].alpha = 2.0;  // Larger α
        outer_data[0].r = 2.0;     // Larger r
        outer_data[0].beta = 2.0;  // Larger f
        outer_data[0].theta = M_PI/2.0;
        outer_data[0].phi_coord = 0.0;
        
        double delta_J = compute_surgical_flux(inner_data, outer_data, true);
        
        // term2 > term1, so ΔJ should be positive
        CHECK(delta_J > 0.0);
        
        // Case 2: α and f decrease outward
        inner_data[0].alpha = 2.0;
        inner_data[0].beta = 2.0;
        
        outer_data[0].alpha = 1.0;
        outer_data[0].beta = 1.0;
        
        delta_J = compute_surgical_flux(inner_data, outer_data, true);
        
        // term2 < term1, so ΔJ should be negative
        CHECK(delta_J < 0.0);
    }
    
    SECTION("Multiple points on surfaces") {
        // Test with multiple points on each surface
        
        int n_points = 4;
        std::vector<MetricData> inner_data(n_points);
        std::vector<MetricData> outer_data(n_points);
        
        double alpha1_avg = 0.0, r1_avg = 0.0, f1_avg = 0.0;
        double alpha2_avg = 0.0, r2_avg = 0.0, f2_avg = 0.0;
        
        for (int i = 0; i < n_points; i++) {
            inner_data[i].alpha = 1.0 + 0.1 * i;
            inner_data[i].r = 1.0;
            inner_data[i].beta = 1.0 + 0.2 * i;  // f
            inner_data[i].theta = M_PI/2.0;
            inner_data[i].phi_coord = 2.0 * M_PI * i / n_points;
            
            outer_data[i].alpha = 2.0 + 0.1 * i;
            outer_data[i].r = 2.0;
            outer_data[i].beta = 2.0 + 0.2 * i;
            outer_data[i].theta = M_PI/2.0;
            outer_data[i].phi_coord = 2.0 * M_PI * i / n_points;
            
            alpha1_avg += inner_data[i].alpha;
            f1_avg += inner_data[i].beta;
            alpha2_avg += outer_data[i].alpha;
            f2_avg += outer_data[i].beta;
        }
        
        alpha1_avg /= n_points;
        f1_avg /= n_points;
        r1_avg = 1.0;
        
        alpha2_avg /= n_points;
        f2_avg /= n_points;
        r2_avg = 2.0;
        
        double delta_J_expected = (1.0/3.0) * 
            (alpha2_avg * r2_avg * r2_avg * f2_avg - 
             alpha1_avg * r1_avg * r1_avg * f1_avg);
        
        double delta_J = compute_surgical_flux(inner_data, outer_data, true);
        
        CHECK_THAT(delta_J, WithinAbs(delta_J_expected, 1e-10));
    }
    
    SECTION("Cumulative angular momentum loss") {
        // Test multiple surgical boundaries
        
        std::vector<std::vector<MetricData>> surgical_boundaries;
        std::vector<bool> axisymmetric_flags;
        
        // Create 3 boundaries
        for (int i = 0; i < 3; i++) {
            std::vector<MetricData> inner(1);
            std::vector<MetricData> outer(1);
            
            inner[0].alpha = 1.0 * (i + 1);
            inner[0].r = 1.0 + i;
            inner[0].beta = 1.0;
            inner[0].theta = M_PI/2.0;
            
            outer[0].alpha = 2.0 * (i + 1);
            outer[0].r = 2.0 + i;
            outer[0].beta = 2.0;
            outer[0].theta = M_PI/2.0;
            
            surgical_boundaries.push_back(inner);
            surgical_boundaries.push_back(outer);
            axisymmetric_flags.push_back(true);
        }
        
        double total_loss = compute_cumulative_angular_momentum_loss(
            surgical_boundaries, axisymmetric_flags);
        
        // Compute expected total
        double total_expected = 0.0;
        for (size_t i = 0; i < surgical_boundaries.size(); i += 2) {
            if (i + 1 < surgical_boundaries.size()) {
                const auto& inner = surgical_boundaries[i];
                const auto& outer = surgical_boundaries[i + 1];
                
                double alpha1 = inner[0].alpha;
                double r1 = inner[0].r;
                double f1 = inner[0].beta;
                
                double alpha2 = outer[0].alpha;
                double r2 = outer[0].r;
                double f2 = outer[0].beta;
                
                total_expected += (1.0/3.0) * (alpha2 * r2 * r2 * f2 - alpha1 * r1 * r1 * f1);
            }
        }
        
        CHECK_THAT(total_loss, WithinAbs(total_expected, 1e-10));
    }
}