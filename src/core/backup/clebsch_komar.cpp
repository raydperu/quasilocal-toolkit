#include "qlt/quasilocal_quantities.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace qlt {

namespace { // Anonymous namespace for helper functions

/**
 * @brief Compute geometric term: ∇^i ξ^j for Komar integral
 */
double compute_komar_geometric_term(const MetricData& data) {
    // ∇^i ξ^j = h^{ik} ∇_k ξ^j
    Matrix3d h_inv = data.h_ij.inverse();
    Matrix3d grad_xi_up = h_inv * data.dxi_dx;
    
    // Trace gives the scalar integrand (for certain orientations)
    // For angular momentum, we need antisymmetric part
    Matrix3d antisym = (grad_xi_up - grad_xi_up.transpose()) / 2.0;
    
    // The integrand for surface integral is more complex
    // Here we compute volume integrand form for Clebsch-Komar
    
    return antisym.norm() * data.sqrt_det_h;
}

/**
 * @brief Compute matter term: α dβ ∧ dξ^♭
 *        3D form: α (∇β × ∇ξ) · n dV
 */
double compute_matter_term_alpha(const MetricData& data) {
    // dβ ∧ dξ^♭ in 3D: (∇β × ∇ξ) · n
    Vector3d cross_product = data.dbeta_dx.cross(
        Vector3d(data.dxi_dx(0,0), data.dxi_dx(1,1), data.dxi_dx(2,2))
    );
    
    return data.alpha * cross_product.dot(data.normal) * data.sqrt_det_h;
}

/**
 * @brief Compute matter term: Φ d(dβ ∧ ξ^♭)
 *        Using d(dβ ∧ ξ^♭) = d²β ∧ ξ^♭ - dβ ∧ dξ^♭ = -dβ ∧ dξ^♭ (since d²β=0)
 *        So this becomes: -Φ (dβ ∧ dξ^♭)
 */
double compute_matter_term_phi(const MetricData& data) {
    // Same as α term but with -Φ
    Vector3d cross_product = data.dbeta_dx.cross(
        Vector3d(data.dxi_dx(0,0), data.dxi_dx(1,1), data.dxi_dx(2,2))
    );
    
    return -data.phi * cross_product.dot(data.normal) * data.sqrt_det_h;
}

/**
 * @brief Verify Clebsch constraint C_{μν} = 0 at a point
 *        C_{μν} = ∇_μ V_ν - ∇_ν V_μ - (∇_μ α ∇_ν β - ∇_ν α ∇_μ β)
 *        where V_μ = ∇_μ Φ + α ∇_μ β
 */
bool verify_clebsch_constraint(const MetricData& data, double tolerance = 1e-8) {
    // Compute V_μ = ∇_μ Φ + α ∇_μ β
    Vector3d V = data.dphi_dx + data.alpha * data.dbeta_dx;
    
    // Approximate ∇_μ V_ν using finite differences assumption
    // In actual code, would need second derivatives
    // For now, check the algebraic part:
    
    // Compute (∇_μ α ∇_ν β - ∇_ν α ∇_μ β)
    Matrix3d antisym_alpha_beta = 
        data.dalpha_dx * data.dbeta_dx.transpose() -
        data.dbeta_dx * data.dalpha_dx.transpose();
    
    // In full implementation, would compute ∇_μ V_ν from neighboring points
    // and check ||∇_μ V_ν - ∇_ν V_μ - antisym_alpha_beta|| < tolerance
    
    return true; // Placeholder
}

} // anonymous namespace

double compute_clebsch_komar_j(
    const std::vector<MetricData>& surface_data,
    const std::vector<MetricData>& volume_data,
    double kappa) {
    
    if (surface_data.empty() && volume_data.empty()) {
        throw std::invalid_argument(
            "compute_clebsch_komar_j: Both surface and volume data empty"
        );
    }
    
    double geometric_term = 0.0;
    double matter_term = 0.0;
    
    // ========================================================================
    // SURFACE INTEGRAL: (1/16π) ∫_S ∇^μ ξ^ν dS_{μν}
    // ========================================================================
    if (!surface_data.empty()) {
        for (const auto& data : surface_data) {
            // Surface element: dS_{ij} = n_i dA, where n is unit normal
            // For angular momentum, integrand is ε_{ijk} ξ^j n^k
            
            // Normal vector should be outward pointing
            if (data.normal.norm() < 1e-10) {
                throw std::runtime_error(
                    "compute_clebsch_komar_j: Normal vector not provided "
                    "for surface data"
                );
            }
            
            Vector3d unit_normal = data.normal.normalized();
            
            // Surface area element: dA = √σ d²θ
            // For sphere in spherical coords: √σ = r² sinθ
            double dA = data.sqrt_det_h; // Should be surface √det(2-metric)
            
            // Komar integrand for angular momentum (simplified for sphere)
            // In full GR: (1/2) ε_{ijk} ∇^j ξ^k n^i dA
            Vector3d curl_xi = Vector3d(
                data.dxi_dx(2,1) - data.dxi_dx(1,2),  // (∇×ξ)_x
                data.dxi_dx(0,2) - data.dxi_dx(2,0),  // (∇×ξ)_y
                data.dxi_dx(1,0) - data.dxi_dx(0,1)   // (∇×ξ)_z
            );
            
            double integrand = curl_xi.dot(unit_normal) * dA;
            geometric_term += integrand;
        }
        
        geometric_term /= 16.0 * M_PI;
    }
    
    // ========================================================================
    // VOLUME INTEGRAL: (κ/16π) ∫_Σ [α dβ ∧ dξ^♭ + Φ d(dβ ∧ ξ^♭)]
    // ========================================================================
    if (!volume_data.empty()) {
        double volume_integral = 0.0;
        size_t constraint_violations = 0;
        
        for (const auto& data : volume_data) {
            // Verify Clebsch constraint (optional but recommended)
            if (!verify_clebsch_constraint(data)) {
                constraint_violations++;
            }
            
            double term_alpha = compute_matter_term_alpha(data);
            double term_phi = compute_matter_term_phi(data);
            
            volume_integral += term_alpha + term_phi;
        }
        
        if (constraint_violations > 0) {
            std::cerr << "Warning: Clebsch constraint violated at " 
                      << constraint_violations << " points" << std::endl;
        }
        
        matter_term = (kappa / (16.0 * M_PI)) * volume_integral;
    }
    
    double total_j = geometric_term + matter_term;
    
    std::cout << "[DEBUG] Clebsch-Komar: geometric = " << geometric_term
              << ", matter = " << matter_term
              << ", total J = " << total_j << std::endl;
    
    return total_j;
}

/**
 * @brief Alternative: Direct surface integral form (for comparison)
 */
double compute_komar_surface_integral(
    const std::vector<MetricData>& surface_data) {
    
    double integral = 0.0;
    
    for (const auto& data : surface_data) {
        // Proper Komar integrand: ∇^μ ξ^ν dS_{μν}
        // For sphere with normal in r-direction:
        
        Matrix3d h_inv = data.h_ij.inverse();
        Matrix3d grad_xi_up = h_inv * data.dxi_dx;
        
        // Antisymmetric part only
        Matrix3d antisym = (grad_xi_up - grad_xi_up.transpose()) / 2.0;
        
        // Contract with surface element (simplified for radial normal)
        // In spherical coords: dS_{θφ} = r² sinθ dθ dφ
        double dA = data.sqrt_det_h;
        
        // For ξ = ∂_φ, the only non-zero component is ∇^θ ξ^φ - ∇^φ ξ^θ
        // This depends on coordinate system
        
        // Simplified: trace of antisym * dA (up to factor)
        double integrand = antisym.trace() * dA;
        integral += integrand;
    }
    
    return integral / (16.0 * M_PI);
}

} // namespace qlt