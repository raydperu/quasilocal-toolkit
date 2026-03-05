#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "qlt/quasilocal_quantities.h"
#include <cmath>
#include <vector>

using namespace qlt;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Axisymmetric formula tests", "[axisymmetric][formula]") {
    SECTION("Basic formula J(r) = (1/3) α r² f") {
        // Test with simple values
        
        double alpha = 2.0;
        double f = 3.0;
        double r = 4.0;
        
        double J_expected = (1.0/3.0) * alpha * r * r * f;
        
        std::vector<double> alpha_vec = {alpha};
        std::vector<double> f_vec = {f};
        std::vector<double> r_vec = {r};
        
        auto J_vec = compute_axisymmetric_j(alpha_vec, f_vec, r_vec);
        
        CHECK(J_vec.size() == 1);
        CHECK_THAT(J_vec[0], WithinAbs(J_expected, 1e-10));
    }
    
    SECTION("Multiple points") {
        int n = 10;
        std::vector<double> alpha_vec(n);
        std::vector<double> f_vec(n);
        std::vector<double> r_vec(n);
        std::vector<double> J_expected_vec(n);
        
        for (int i = 0; i < n; i++) {
            alpha_vec[i] = 0.1 * (i + 1);
            f_vec[i] = 0.2 * (i + 1);
            r_vec[i] = 1.0 + i;
            J_expected_vec[i] = (1.0/3.0) * alpha_vec[i] * r_vec[i] * r_vec[i] * f_vec[i];
        }
        
        auto J_computed_vec = compute_axisymmetric_j(alpha_vec, f_vec, r_vec);
        
        CHECK(J_computed_vec.size() == n);
        
        for (int i = 0; i < n; i++) {
            CHECK_THAT(J_computed_vec[i], WithinAbs(J_expected_vec[i], 1e-10));
        }
    }
    
    SECTION("Zero cases") {
        // Test cases where J should be zero
        
        SECTION("α = 0") {
            std::vector<double> alpha_vec = {0.0, 0.0, 0.0};
            std::vector<double> f_vec = {1.0, 2.0, 3.0};
            std::vector<double> r_vec = {1.0, 2.0, 3.0};
            
            auto J_vec = compute_axisymmetric_j(alpha_vec, f_vec, r_vec);
            
            for (double J : J_vec) {
                CHECK_THAT(J, WithinAbs(0.0, 1e-10));
            }
        }
        
        SECTION("f = 0") {
            std::vector<double> alpha_vec = {1.0, 2.0, 3.0};
            std::vector<double> f_vec = {0.0, 0.0, 0.0};
            std::vector<double> r_vec = {1.0, 2.0, 3.0};
            
            auto J_vec = compute_axisymmetric_j(alpha_vec, f_vec, r_vec);
            
            for (double J : J_vec) {
                CHECK_THAT(J, WithinAbs(0.0, 1e-10));
            }
        }
        
        SECTION("r = 0") {
            std::vector<double> alpha_vec = {1.0};
            std::vector<double> f_vec = {2.0};
            std::vector<double> r_vec = {0.0};
            
            auto J_vec = compute_axisymmetric_j(alpha_vec, f_vec, r_vec);
            
            CHECK(J_vec.size() == 1);
            CHECK_THAT(J_vec[0], WithinAbs(0.0, 1e-10));
        }
    }
    
    SECTION("Negative values") {
        // Test with negative α, f, r
        
        SECTION("Negative α") {
            std::vector<double> alpha_vec = {-1.0, -2.0, -3.0};
            std::vector<double> f_vec = {1.0, 1.0, 1.0};
            std::vector<double> r_vec = {1.0, 2.0, 3.0};
            
            auto J_vec = compute_axisymmetric_j(alpha_vec, f_vec, r_vec);
            
            for (int i = 0; i < 3; i++) {
                double J_expected = (1.0/3.0) * alpha_vec[i] * r_vec[i] * r_vec[i] * f_vec[i];
                CHECK_THAT(J_vec[i], WithinAbs(J_expected, 1e-10));
                CHECK(J_vec[i] < 0.0);  // Should be negative
            }
        }
        
        SECTION("Negative f") {
            std::vector<double> alpha_vec = {1.0, 1.0, 1.0};
            std::vector<double> f_vec = {-1.0, -2.0, -3.0};
            std::vector<double> r_vec = {1.0, 2.0, 3.0};
            
            auto J_vec = compute_axisymmetric_j(alpha_vec, f_vec, r_vec);
            
            for (int i = 0; i < 3; i++) {
                double J_expected = (1.0/3.0) * alpha_vec[i] * r_vec[i] * r_vec[i] * f_vec[i];
                CHECK_THAT(J_vec[i], WithinAbs(J_expected, 1e-10));
                CHECK(J_vec[i] < 0.0);
            }
        }
        
        SECTION("Negative r (should throw)") {
            std::vector<double> alpha_vec = {1.0};
            std::vector<double> f_vec = {1.0};
            std::vector<double> r_vec = {-1.0};  // Negative radius
            
            CHECK_THROWS_AS(
                compute_axisymmetric_j(alpha_vec, f_vec, r_vec),
                std::domain_error
            );
        }
    }
    
    SECTION("Large values") {
        // Test with very large/small numbers
        
        SECTION("Large radius") {
            double alpha = 1.0;
            double f = 1.0;
            double r = 1.0e6;  // Large radius
            
            std::vector<double> alpha_vec = {alpha};
            std::vector<double> f_vec = {f};
            std::vector<double> r_vec = {r};
            
            auto J_vec = compute_axisymmetric_j(alpha_vec, f_vec, r_vec);
            
            double J_expected = (1.0/3.0) * alpha * r * r * f;
            CHECK_THAT(J_vec[0], WithinRel(J_expected, 1e-10));
        }
        
        SECTION("Small values") {
            double alpha = 1.0e-6;
            double f = 1.0e-6;
            double r = 1.0e-3;
            
            std::vector<double> alpha_vec = {alpha};
            std::vector<double> f_vec = {f};
            std::vector<double> r_vec = {r};
            
            auto J_vec = compute_axisymmetric_j(alpha_vec, f_vec, r_vec);
            
            double J_expected = (1.0/3.0) * alpha * r * r * f;
            CHECK_THAT(J_vec[0], WithinRel(J_expected, 1e-10));
        }
    }
    
    SECTION("Array size mismatch") {
        // Test error handling for mismatched array sizes
        
        std::vector<double> alpha_vec = {1.0, 2.0};
        std::vector<double> f_vec = {1.0};  // Different size
        std::vector<double> r_vec = {1.0, 2.0};
        
        CHECK_THROWS_AS(
            compute_axisymmetric_j(alpha_vec, f_vec, r_vec),
            std::invalid_argument
        );
    }
    
    SECTION("Empty arrays") {
        // Test with empty input
        
        std::vector<double> alpha_vec = {};
        std::vector<double> f_vec = {};
        std::vector<double> r_vec = {};
        
        auto J_vec = compute_axisymmetric_j(alpha_vec, f_vec, r_vec);
        
        CHECK(J_vec.empty());
    }
}