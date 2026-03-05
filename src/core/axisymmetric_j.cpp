#include "quasilocal_quantities.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace qlt {

std::vector<double> compute_axisymmetric_j(
    const std::vector<double>& alpha,
    const std::vector<double>& f,
    const std::vector<double>& r) {
    
    // Validate input sizes
    if (alpha.size() != f.size() || f.size() != r.size()) {
        throw std::invalid_argument(
            "compute_axisymmetric_j: All input vectors must have same size. "
            "Got alpha=" + std::to_string(alpha.size()) +
            ", f=" + std::to_string(f.size()) +
            ", r=" + std::to_string(r.size())
        );
    }
    
    // Check for empty input
    if (alpha.empty()) {
        return {};
    }
    
    // Apply Theorem 4.3: J(r) = (1/3) α(r) r² f(r)
    std::vector<double> J(alpha.size());
    for (size_t i = 0; i < alpha.size(); ++i) {
        // Validate inputs
        if (r[i] < 0.0) {
            throw std::domain_error(
                "compute_axisymmetric_j: Radial coordinate r[" + 
                std::to_string(i) + "] = " + std::to_string(r[i]) + 
                " is negative"
            );
        }
        
        // Main formula
        J[i] = (1.0 / 3.0) * alpha[i] * r[i] * r[i] * f[i];
        
        // Debug output for first few points
        if (i < 3) {
            std::cout << "[DEBUG] Axisymmetric J: r=" << r[i] 
                      << ", α=" << alpha[i] 
                      << ", f=" << f[i] 
                      << ", J=" << J[i] << std::endl;
        }
    }
    
    return J;
}

/**
 * @brief Axisymmetric surgical flux (Corollary 4.6)
 *        ΔJ = (1/3)[α(r₂)r₂²f(r₂) - α(r₁)r₁²f(r₁)]
 */
double compute_axisymmetric_surgical_flux(
    double alpha1, double r1, double f1,
    double alpha2, double r2, double f2) {
    
    if (r1 < 0.0 || r2 < 0.0) {
        throw std::domain_error(
            "compute_axisymmetric_surgical_flux: Negative radii"
        );
    }
    
    double term1 = alpha1 * r1 * r1 * f1;
    double term2 = alpha2 * r2 * r2 * f2;
    
    return (1.0 / 3.0) * (term2 - term1);
}

/**
 * @brief Extract α(r) from axisymmetric data assuming β = φ + f(r)cosθ
 *        and data is given on θ = π/2 slice (equatorial plane)
 */
std::vector<double> extract_alpha_from_surface(
    const std::vector<MetricData>& surface_data) {
    
    std::vector<double> alpha_vals;
    alpha_vals.reserve(surface_data.size());
    
    for (const auto& data : surface_data) {
        // Check if near equatorial plane (for axisymmetry)
        if (std::abs(data.theta - M_PI/2.0) > 1e-6) {
            std::cerr << "Warning: Data not at θ=π/2. Using anyway." << std::endl;
        }
        alpha_vals.push_back(data.alpha);
    }
    
    return alpha_vals;
}

/**
 * @brief Extract f(r) = β - φ from axisymmetric data
 *        Assumes β = φ + f(r)cosθ
 */
std::vector<double> extract_f_from_surface(
    const std::vector<MetricData>& surface_data) {
    
    std::vector<double> f_vals;
    f_vals.reserve(surface_data.size());
    
    for (const auto& data : surface_data) {
        // For β = φ + f(r)cosθ, we have f(r) = (β - φ)/cosθ
        double cos_theta = std::cos(data.theta);
        
        if (std::abs(cos_theta) < 1e-10) {
            // At poles, f(r) is ill-defined from this formula
            // Use average of neighbors or special handling
            f_vals.push_back(0.0);
        } else {
            f_vals.push_back((data.beta - data.phi_coord) / cos_theta);
        }
    }
    
    return f_vals;
}

/**
 * @brief Extract radial coordinates from surface data
 */
std::vector<double> extract_r_from_surface(
    const std::vector<MetricData>& surface_data) {
    
    std::vector<double> r_vals;
    r_vals.reserve(surface_data.size());
    
    for (const auto& data : surface_data) {
        r_vals.push_back(data.r);
    }
    
    return r_vals;
}

} // namespace qlt
