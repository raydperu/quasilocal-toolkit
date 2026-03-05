#include "quasilocal_quantities.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace qlt {

namespace { // Helper functions

/**
 * @brief Compute extrinsic curvature of 2-surface embedded in 3-slice
 *        k = h^{ij} ∇_i n_j (trace of extrinsic curvature)
 */
double compute_surface_extrinsic_curvature(const MetricData& data) {
    // For a 2-surface with normal n, the extrinsic curvature is
    // k = ∇_i n^i = h^{ij} ∇_i n_j
    
    if (data.normal.norm() < 1e-10) {
        throw std::runtime_error(
            "compute_surface_extrinsic_curvature: Normal vector not provided"
        );
    }
    
    Vector3d unit_normal = data.normal.normalized();
    
    // In full GR, we'd need Christoffel symbols to compute ∇_i n_j
    // For now, approximate using the metric and normal derivatives
    // k ≈ h^{ij} (∂_i n_j - Γ^l_{ij} n_l)
    
    // Simplified: Assume normal is radial in spherical coordinates
    // For sphere of radius r: k = 2/r
    
    if (data.r > 1e-10) {
        // Spherical approximation
        return 2.0 / data.r;
    } else {
        return 0.0;
    }
}

/**
 * @brief Compute reference (flat space) extrinsic curvature
 *        For sphere in flat space: k₀ = 2/r
 */
double compute_reference_extrinsic_curvature(const MetricData& data) {
    if (data.r > 1e-10) {
        return 2.0 / data.r;
    } else {
        return 0.0;
    }
}

/**
 * @brief Compute surface area element √σ
 *        For sphere: √σ = r² sinθ
 */
double compute_surface_area_element(const MetricData& data) {
    // sqrt_det_h should be the 2-metric determinant for surface
    // If data.sqrt_det_h is already surface element, use it
    // Otherwise compute from metric
    
    if (data.sqrt_det_h > 0) {
        return data.sqrt_det_h;
    } else {
        // Fallback: spherical coordinates
        if (data.r > 0) {
            return data.r * data.r * std::sin(data.theta);
        } else {
            return 1.0;
        }
    }
}

} // anonymous namespace

double compute_quasilocal_mass(
    const std::vector<MetricData>& surface_data,
    bool use_reference) {
    
    if (surface_data.empty()) {
        throw std::invalid_argument(
            "compute_quasilocal_mass: Surface data is empty"
        );
    }
    
    double total_mass = 0.0;
    double total_area = 0.0;
    
    for (const auto& data : surface_data) {
        // Compute extrinsic curvature of surface
        double k = compute_surface_extrinsic_curvature(data);
        
        // Compute reference curvature (flat space)
        double k0 = use_reference ? compute_reference_extrinsic_curvature(data) : 0.0;
        
        // Surface area element
        double dA = compute_surface_area_element(data);
        
        // Brown-York style quasilocal mass
        // M_S = -(1/8π) ∫_S (k - k₀) √σ d²θ
        double mass_contribution = -(k - k0) * dA;
        
        total_mass += mass_contribution;
        total_area += dA;
        
        // Debug output
        if (&data == &surface_data[0]) { // First point only
            std::cout << "[DEBUG] Quasilocal mass: k = " << k
                      << ", k0 = " << k0
                      << ", dA = " << dA
                      << ", contribution = " << mass_contribution << std::endl;
        }
    }
    
    total_mass /= (8.0 * M_PI);
    
    std::cout << "[DEBUG] Total quasilocal mass = " << total_mass
              << ", total area = " << total_area << std::endl;
    
    return total_mass;
}

/**
 * @brief Compute ADM mass at infinity (for comparison)
 *        For asymptotically flat spacetimes
 */
double compute_adm_mass_approximation(
    const std::vector<MetricData>& surface_data) {
    
    // Simplified: Assuming spherical surface at large r
    // ADM mass ≈ (r/2) (1 - 1/√g_rr) for Schwarzschild
    
    if (surface_data.empty()) {
        return 0.0;
    }
    
    // Use the data point with largest r
    double max_r = 0.0;
    double grr_at_max_r = 1.0;
    
    for (const auto& data : surface_data) {
        if (data.r > max_r) {
            max_r = data.r;
            // g_rr component (assuming spherical coordinates)
            grr_at_max_r = data.h_ij(0, 0); // First diagonal component
        }
    }
    
    if (max_r < 1e-10) {
        return 0.0;
    }
    
    // For Schwarzschild: g_rr = 1/(1 - 2M/r)
    // So M ≈ (r/2) (1 - 1/g_rr)
    double adm_mass = (max_r / 2.0) * (1.0 - 1.0 / grr_at_max_r);
    
    return std::max(0.0, adm_mass); // Mass should be non-negative
}

/**
 * @brief Compute mass aspect function M(r) for spherical symmetry
 */
std::vector<double> compute_mass_aspect(
    const std::vector<MetricData>& surface_data_at_radii) {
    
    std::vector<double> mass_aspect;
    mass_aspect.reserve(surface_data_at_radii.size());
    
    for (const auto& data : surface_data_at_radii) {
        // For spherical symmetry: M(r) = (r/2)(1 - h^{rr})
        // where h^{rr} is inverse radial metric component
        
        Matrix3d h_inv = data.h_ij.inverse();
        double h_rr_inv = h_inv(0, 0); // Assuming (r,θ,φ) ordering
        
        double mass = (data.r / 2.0) * (1.0 - 1.0 / h_rr_inv);
        mass_aspect.push_back(std::max(0.0, mass));
    }
    
    return mass_aspect;
}

} // namespace qlt
