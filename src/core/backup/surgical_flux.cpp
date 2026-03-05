#include "qlt/quasilocal_quantities.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace qlt {

namespace { // Helper functions

/**
 * @brief Compute integrand for general surgical flux
 *        Theorem 4.4: α dβ ∧ ξ^♭
 */
double compute_surgical_flux_integrand(const MetricData& data) {
    // α dβ ∧ ξ^♭ in 3D form
    // For surface integral: α (∇β × ξ) · n dA
    
    if (data.normal.norm() < 1e-10) {
        throw std::runtime_error(
            "compute_surgical_flux_integrand: Normal vector not provided"
        );
    }
    
    Vector3d unit_normal = data.normal.normalized();
    
    // ξ^♭ is the 1-form dual to ξ: ξ_i = g_{ij} ξ^j
    Vector3d xi_flat = data.h_ij * data.xi;
    
    // Cross product: ∇β × ξ
    Vector3d cross_product = data.dbeta_dx.cross(xi_flat);
    
    // Integrand: α (∇β × ξ) · n
    double integrand = data.alpha * cross_product.dot(unit_normal);
    
    // Multiply by surface area element
    double dA = data.sqrt_det_h; // Should be √σ for surface
    
    return integrand * dA;
}

} // anonymous namespace

double compute_surgical_flux(
    const std::vector<MetricData>& inner_data,
    const std::vector<MetricData>& outer_data,
    bool axisymmetric) {
    
    if (inner_data.empty() && outer_data.empty()) {
        return 0.0;
    }
    
    // ========================================================================
    // AXISYMMETRIC CASE (Corollary 4.6)
    // ΔJ = (1/3)[α(r₂)r₂²f(r₂) - α(r₁)r₁²f(r₁)]
    // ========================================================================
    if (axisymmetric) {
        if (inner_data.empty() || outer_data.empty()) {
            throw std::invalid_argument(
                "compute_surgical_flux: Axisymmetric case requires "
                "both inner and outer surface data"
            );
        }
        
        // Extract average values from surfaces
        double alpha1_avg = 0.0, r1_avg = 0.0, f1_avg = 0.0;
        double alpha2_avg = 0.0, r2_avg = 0.0, f2_avg = 0.0;
        
        for (const auto& data : inner_data) {
            alpha1_avg += data.alpha;
            r1_avg += data.r;
            // Extract f(r) = (β - φ)/cosθ
            double cos_theta = std::cos(data.theta);
            if (std::abs(cos_theta) > 1e-10) {
                f1_avg += (data.beta - data.phi_coord) / cos_theta;
            }
        }
        
        for (const auto& data : outer_data) {
            alpha2_avg += data.alpha;
            r2_avg += data.r;
            double cos_theta = std::cos(data.theta);
            if (std::abs(cos_theta) > 1e-10) {
                f2_avg += (data.beta - data.phi_coord) / cos_theta;
            }
        }
        
        size_t n1 = inner_data.size();
        size_t n2 = outer_data.size();
        
        alpha1_avg /= n1; r1_avg /= n1; f1_avg /= n1;
        alpha2_avg /= n2; r2_avg /= n2; f2_avg /= n2;
        
        // Apply Corollary 4.6
        double delta_j = (1.0/3.0) * 
            (alpha2_avg * r2_avg * r2_avg * f2_avg - 
             alpha1_avg * r1_avg * r1_avg * f1_avg);
        
        std::cout << "[DEBUG] Axisymmetric surgical flux:" << std::endl;
        std::cout << "  Inner: α=" << alpha1_avg << ", r=" << r1_avg 
                  << ", f=" << f1_avg << std::endl;
        std::cout << "  Outer: α=" << alpha2_avg << ", r=" << r2_avg 
                  << ", f=" << f2_avg << std::endl;
        std::cout << "  ΔJ = " << delta_j << std::endl;
        
        return delta_j;
    }
    
    // ========================================================================
    // GENERAL CASE (Theorem 4.4)
    // ΔJ = -(1/8π) ∫_∂V α dβ ∧ ξ^♭
    // ========================================================================
    else {
        double flux_integral = 0.0;
        
        // Outer boundary contributes positively (outward normal)
        for (const auto& data : outer_data) {
            flux_integral += compute_surgical_flux_integrand(data);
        }
        
        // Inner boundary contributes negatively (inward normal, but
        // normal should be outward from excised region, so flip sign)
        for (const auto& data : inner_data) {
            // For inner surface, normal points into excised region
            // So contribution is negative
            flux_integral -= compute_surgical_flux_integrand(data);
        }
        
        double delta_j = -flux_integral / (8.0 * M_PI);
        
        std::cout << "[DEBUG] General surgical flux:" << std::endl;
        std::cout << "  Flux integral = " << flux_integral << std::endl;
        std::cout << "  ΔJ = " << delta_j << std::endl;
        
        return delta_j;
    }
}

/**
 * @brief Compute cumulative angular momentum loss through multiple surgeries
 */
double compute_cumulative_angular_momentum_loss(
    const std::vector<std::vector<MetricData>>& surgical_boundaries,
    const std::vector<bool>& axisymmetric_flags) {
    
    if (surgical_boundaries.size() < 2) {
        throw std::invalid_argument(
            "compute_cumulative_angular_momentum_loss: "
            "Need at least 2 boundaries"
        );
    }
    
    double total_loss = 0.0;
    
    // Compute flux between consecutive boundaries
    for (size_t i = 0; i < surgical_boundaries.size() - 1; i += 2) {
        // Assume boundaries come in pairs: (inner, outer)
        if (i + 1 >= surgical_boundaries.size()) {
            break;
        }
        
        bool axisym = (i < axisymmetric_flags.size()) ? 
                     axisymmetric_flags[i] : false;
        
        double delta_j = compute_surgical_flux(
            surgical_boundaries[i],    // inner
            surgical_boundaries[i+1],  // outer
            axisym
        );
        
        total_loss += delta_j;
        
        std::cout << "[DEBUG] Surgery " << i/2 << ": ΔJ = " << delta_j
                  << ", cumulative = " << total_loss << std::endl;
    }
    
    return total_loss;
}

/**
 * @brief Verify angular momentum conservation across surgery
 *        Check: J_after = J_before + ΔJ (should be conserved)
 */
bool verify_surgery_conservation(
    double j_before,           // J(λ⁻)
    double j_after,            // J(λ⁺)
    const std::vector<MetricData>& inner_boundary,
    const std::vector<MetricData>& outer_boundary,
    double tolerance = 1e-6) {
    
    double delta_j_computed = compute_surgical_flux(
        inner_boundary, outer_boundary, false
    );
    
    double j_predicted = j_before + delta_j_computed;
    double error = std::abs(j_after - j_predicted);
    
    std::cout << "[VERIFICATION] Surgery conservation check:" << std::endl;
    std::cout << "  J_before = " << j_before << std::endl;
    std::cout << "  J_after = " << j_after << std::endl;
    std::cout << "  ΔJ_computed = " << delta_j_computed << std::endl;
    std::cout << "  J_predicted = " << j_predicted << std::endl;
    std::cout << "  Error = " << error << " (tolerance = " << tolerance << ")" << std::endl;
    
    return error < tolerance;
}

} // namespace qlt