// surface_integral.cpp 
#include "surface_integral.h" 
#include "coordinate_systems.h" 
 
namespace qlt { 
namespace grid { 
 
#include "surface_integral.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>

namespace qlt {
namespace grid {

namespace { // Helper functions

double compute_midpoint_integral(const std::vector<double>& values,
                                 const std::vector<double>& weights) {
    if (values.size() != weights.size()) {
        throw std::invalid_argument("Values and weights must have same size");
    }
    
    double integral = 0.0;
    for (size_t i = 0; i < values.size(); i++) {
        integral += values[i] * weights[i];
    }
    return integral;
}

double compute_trapezoidal_integral(const std::vector<double>& values,
                                    const std::vector<double>& weights) {
    // Simplified trapezoidal rule
    return compute_midpoint_integral(values, weights);
}

std::vector<double> compute_connectivity_based_weights(
    const std::vector<MetricData>& surface_data,
    const std::vector<std::array<size_t, 3>>& faces) {
    
    std::vector<double> weights(surface_data.size(), 0.0);
    
    for (const auto& face : faces) {
        // Get vertices
        const auto& v0 = surface_data[face[0]];
        const auto& v1 = surface_data[face[1]];
        const auto& v2 = surface_data[face[2]];
        
        // Compute area of triangle
        // Convert to Cartesian for area calculation
        auto pos0 = CoordinateSystems::spherical_to_cartesian(v0.r, v0.theta, v0.phi_coord);
        auto pos1 = CoordinateSystems::spherical_to_cartesian(v1.r, v1.theta, v1.phi_coord);
        auto pos2 = CoordinateSystems::spherical_to_cartesian(v2.r, v2.theta, v2.phi_coord);
        
        Vector3d p0(pos0[0], pos0[1], pos0[2]);
        Vector3d p1(pos1[0], pos1[1], pos1[2]);
        Vector3d p2(pos2[0], pos2[1], pos2[2]);
        
        Vector3d v = p1 - p0;
        Vector3d w = p2 - p0;
        double area = 0.5 * v.cross(w).norm();
        
        // Distribute area to vertices (equal share)
        weights[face[0]] += area / 3.0;
        weights[face[1]] += area / 3.0;
        weights[face[2]] += area / 3.0;
    }
    
    return weights;
}

} // anonymous namespace

double SurfaceIntegral::integrate_scalar(
    const std::vector<MetricData>& surface_data,
    const std::vector<double>& scalar_values,
    IntegrationMethod method) {
    
    if (surface_data.size() != scalar_values.size()) {
        throw std::invalid_argument(
            "Surface data and scalar values must have same size");
    }
    
    // Compute area elements
    std::vector<double> area_elements = compute_area_elements(surface_data);
    
    switch (method) {
        case IntegrationMethod::MIDPOINT_RULE:
            return compute_midpoint_integral(scalar_values, area_elements);
        
        case IntegrationMethod::TRAPEZOIDAL_RULE:
            return compute_trapezoidal_integral(scalar_values, area_elements);
        
        case IntegrationMethod::SPHERICAL_QUADRATURE: {
            // Special quadrature for spheres
            double integral = 0.0;
            for (size_t i = 0; i < surface_data.size(); i++) {
                double dA = area_elements[i];
                // For sphere, weight by sinθ for proper quadrature
                double weight = sin(surface_data[i].theta);
                integral += scalar_values[i] * dA * weight;
            }
            return integral;
        }
        
        default:
            throw std::runtime_error("Unsupported integration method");
    }
}

double SurfaceIntegral::integrate_vector_flux(
    const std::vector<MetricData>& surface_data,
    const std::vector<Vector3d>& vector_field,
    IntegrationMethod method) {
    
    if (surface_data.size() != vector_field.size()) {
        throw std::invalid_argument(
            "Surface data and vector field must have same size");
    }
    
    // Compute dot product with normal at each point
    std::vector<double> normal_flux(surface_data.size());
    for (size_t i = 0; i < surface_data.size(); i++) {
        if (surface_data[i].normal.norm() < 1e-10) {
            throw std::runtime_error(
                "Normal vector not provided for flux calculation");
        }
        
        Vector3d unit_normal = surface_data[i].normal.normalized();
        normal_flux[i] = vector_field[i].dot(unit_normal);
    }
    
    // Integrate scalar flux
    return integrate_scalar(surface_data, normal_flux, method);
}

double SurfaceIntegral::compute_komar_surface_integral(
    const std::vector<MetricData>& surface_data) {
    
    double integral = 0.0;
    
    for (const auto& data : surface_data) {
        // Komar integrand: ∇^μ ξ^ν dS_{μν}
        // Simplified for spherical surface
        
        if (data.normal.norm() < 1e-10) {
            throw std::runtime_error("Normal vector not provided");
        }
        
        // Compute curl of ξ (∇ × ξ)
        Vector3d curl_xi(
            data.dxi_dx(2,1) - data.dxi_dx(1,2),  // (∇×ξ)_x
            data.dxi_dx(0,2) - data.dxi_dx(2,0),  // (∇×ξ)_y
            data.dxi_dx(1,0) - data.dxi_dx(0,1)   // (∇×ξ)_z
        );
        
        Vector3d unit_normal = data.normal.normalized();
        double integrand = curl_xi.dot(unit_normal);
        
        // Area element
        double dA = data.sqrt_det_h;
        
        integral += integrand * dA;
    }
    
    return integral / (16.0 * M_PI);
}

double SurfaceIntegral::compute_surface_area(
    const std::vector<MetricData>& surface_data) {
    
    std::vector<double> area_elements = compute_area_elements(surface_data);
    return std::accumulate(area_elements.begin(), area_elements.end(), 0.0);
}

std::array<double, 3> SurfaceIntegral::compute_surface_centroid(
    const std::vector<MetricData>& surface_data) {
    
    std::array<double, 3> centroid = {0.0, 0.0, 0.0};
    double total_area = 0.0;
    
    for (const auto& data : surface_data) {
        // Convert to Cartesian
        auto cartesian = CoordinateSystems::spherical_to_cartesian(
            data.r, data.theta, data.phi_coord);
        
        double dA = data.sqrt_det_h;
        
        centroid[0] += cartesian[0] * dA;
        centroid[1] += cartesian[1] * dA;
        centroid[2] += cartesian[2] * dA;
        total_area += dA;
    }
    
    if (total_area > 1e-10) {
        centroid[0] /= total_area;
        centroid[1] /= total_area;
        centroid[2] /= total_area;
    }
    
    return centroid;
}

std::vector<Vector3d> SurfaceIntegral::compute_surface_normals(
    const std::vector<MetricData>& surface_data) {
    
    // Generate connectivity if needed
    auto faces = generate_face_connectivity(surface_data);
    
    std::vector<Vector3d> normals(surface_data.size(), Vector3d::Zero());
    std::vector<int> face_count(surface_data.size(), 0);
    
    for (const auto& face : faces) {
        // Get vertices
        const auto& v0 = surface_data[face[0]];
        const auto& v1 = surface_data[face[1]];
        const auto& v2 = surface_data[face[2]];
        
        // Convert to Cartesian
        auto p0 = CoordinateSystems::spherical_to_cartesian(v0.r, v0.theta, v0.phi_coord);
        auto p1 = CoordinateSystems::spherical_to_cartesian(v1.r, v1.theta, v1.phi_coord);
        auto p2 = CoordinateSystems::spherical_to_cartesian(v2.r, v2.theta, v2.phi_coord);
        
        Vector3d vec0(p0[0], p0[1], p0[2]);
        Vector3d vec1(p1[0], p1[1], p1[2]);
        Vector3d vec2(p2[0], p2[1], p2[2]);
        
        // Compute face normal
        Vector3d v = vec1 - vec0;
        Vector3d w = vec2 - vec0;
        Vector3d face_normal = v.cross(w).normalized();
        
        // Add to vertex normals
        normals[face[0]] += face_normal;
        normals[face[1]] += face_normal;
        normals[face[2]] += face_normal;
        
        face_count[face[0]]++;
        face_count[face[1]]++;
        face_count[face[2]]++;
    }
    
    // Average and normalize
    for (size_t i = 0; i < normals.size(); i++) {
        if (face_count[i] > 0) {
            normals[i] /= face_count[i];
            normals[i].normalize();
        } else {
            // If no faces, use radial direction
            auto cartesian = CoordinateSystems::spherical_to_cartesian(
                surface_data[i].r, surface_data[i].theta, surface_data[i].phi_coord);
            normals[i] = Vector3d(cartesian[0], cartesian[1], cartesian[2]).normalized();
        }
    }
    
    return normals;
}

std::vector<double> SurfaceIntegral::compute_surface_curvature(
    const std::vector<MetricData>& surface_data) {
    
    std::vector<double> curvature(surface_data.size(), 0.0);
    
    // For sphere, curvature = 1/r²
    for (size_t i = 0; i < surface_data.size(); i++) {
        if (surface_data[i].r > 1e-10) {
            curvature[i] = 1.0 / (surface_data[i].r * surface_data[i].r);
        }
    }
    
    return curvature;
}

bool SurfaceIntegral::is_closed_surface(
    const std::vector<MetricData>& surface_data) {
    
    // Simple check: if normals are consistently outward/inward
    // For sphere, this is true
    
    if (surface_data.empty()) return false;
    
    // Check if all points have similar normal directions
    Vector3d first_normal;
    if (surface_data[0].normal.norm() > 1e-10) {
        first_normal = surface_data[0].normal.normalized();
    } else {
        // Compute normal
        auto normals = compute_surface_normals(surface_data);
        first_normal = normals[0];
    }
    
    for (size_t i = 1; i < surface_data.size(); i++) {
        Vector3d current_normal;
        if (surface_data[i].normal.norm() > 1e-10) {
            current_normal = surface_data[i].normal.normalized();
        } else {
            auto normals = compute_surface_normals(surface_data);
            current_normal = normals[i];
        }
        
        // Check if normals point in similar direction
        double dot_product = first_normal.dot(current_normal);
        if (dot_product < 0.5) {  // Arbitrary threshold
            return false;
        }
    }
    
    return true;
}

std::vector<MetricData> SurfaceIntegral::extract_spherical_surface(
    const std::vector<std::vector<std::vector<MetricData>>>& grid_data,
    double radius,
    size_t num_theta,
    size_t num_phi) {
    
    // Simplified extraction - would use interpolation in real implementation
    std::vector<MetricData> surface_data;
    
    // Generate spherical coordinates
    auto spherical_grid = CoordinateSystems::generate_spherical_grid(
        radius, radius, 1, num_theta, num_phi);
    
    // Convert to Cartesian and find nearest grid points
    for (size_t i = 0; i < num_theta; i++) {
        for (size_t j = 0; j < num_phi; j++) {
            auto spherical_coord = spherical_grid[0][i][j];
            auto cartesian = CoordinateSystems::spherical_to_cartesian(
                spherical_coord[0], spherical_coord[1], spherical_coord[2]);
            
            // Find nearest grid point (simplified)
            // In production, would use trilinear interpolation
            
            MetricData data;
            data.r = radius;
            data.theta = spherical_coord[1];
            data.phi_coord = spherical_coord[2];
            
            // Set normal (radial outward)
            data.normal = Vector3d(
                sin(data.theta) * cos(data.phi_coord),
                sin(data.theta) * sin(data.phi_coord),
                cos(data.theta)
            );
            
            surface_data.push_back(data);
        }
    }
    
    return surface_data;
}

std::vector<MetricData> SurfaceIntegral::refine_surface(
    const std::vector<MetricData>& surface_data,
    int refinement_level) {
    
    if (refinement_level <= 0) {
        return surface_data;
    }
    
    // Simple refinement: add midpoints
    std::vector<MetricData> refined_data;
    
    // For now, just return original (would implement subdivision in production)
    refined_data = surface_data;
    
    // Apply smoothing if requested
    if (refinement_level > 1) {
        refined_data = smooth_surface(refined_data, refinement_level - 1);
    }
    
    return refined_data;
}

std::vector<MetricData> SurfaceIntegral::smooth_surface(
    const std::vector<MetricData>& surface_data,
    int smoothing_iterations,
    double smoothing_factor) {
    
    if (smoothing_iterations <= 0) {
        return surface_data;
    }
    
    std::vector<MetricData> smoothed_data = surface_data;
    
    // Simple Laplacian smoothing
    for (int iter = 0; iter < smoothing_iterations; iter++) {
        std::vector<MetricData> temp_data = smoothed_data;
        
        // For each point, average with neighbors
        for (size_t i = 0; i < smoothed_data.size(); i++) {
            // Simplified: just adjust radius slightly
            // In production, would use actual mesh connectivity
            
            double radius_change = 0.0;
            int neighbor_count = 0;
            
            // Check nearby points (simplified)
            for (size_t j = 0; j < smoothed_data.size(); j++) {
                if (i != j) {
                    double dr = smoothed_data[j].r - smoothed_data[i].r;
                    double angular_dist = std::acos(
                        sin(smoothed_data[i].theta) * sin(smoothed_data[j].theta) *
                        cos(smoothed_data[i].phi_coord - smoothed_data[j].phi_coord) +
                        cos(smoothed_data[i].theta) * cos(smoothed_data[j].theta)
                    );
                    
                    if (angular_dist < 0.3) {  // Arbitrary threshold
                        radius_change += dr;
                        neighbor_count++;
                    }
                }
            }
            
            if (neighbor_count > 0) {
                temp_data[i].r += smoothing_factor * (radius_change / neighbor_count);
            }
        }
        
        smoothed_data = temp_data;
    }
    
    return smoothed_data;
}

double SurfaceIntegral::compute_surgical_flux_integral(
    const std::vector<MetricData>& surface_data) {
    
    double integral = 0.0;
    
    for (const auto& data : surface_data) {
        // Surgical flux integrand: α (∇β × ξ) · n
        if (data.normal.norm() < 1e-10) {
            throw std::runtime_error("Normal vector not provided");
        }
        
        Vector3d unit_normal = data.normal.normalized();
        Vector3d xi_flat = data.h_ij * data.xi;
        
        Vector3d cross_product = data.dbeta_dx.cross(xi_flat);
        double integrand = data.alpha * cross_product.dot(unit_normal);
        
        double dA = data.sqrt_det_h;
        integral += integrand * dA;
    }
    
    return integral;
}

// Private helper methods
std::vector<double> SurfaceIntegral::compute_area_elements(
    const std::vector<MetricData>& surface_data) {
    
    std::vector<double> area_elements(surface_data.size());
    
    for (size_t i = 0; i < surface_data.size(); i++) {
        if (surface_data[i].sqrt_det_h > 0.0) {
            // Use provided area element
            area_elements[i] = surface_data[i].sqrt_det_h;
        } else {
            // Compute from metric
            double det_h = surface_data[i].h_ij.determinant();
            area_elements[i] = std::sqrt(std::max(0.0, det_h));
        }
    }
    
    return area_elements;
}

std::vector<std::array<size_t, 3>> SurfaceIntegral::generate_face_connectivity(
    const std::vector<MetricData>& surface_data) {
    
    std::vector<std::array<size_t, 3>> faces;
    
    // Simplified: assume points are arranged in θ-φ grid
    // Find dimensions
    size_t n_theta = 0, n_phi = 0;
    
    // Count unique θ values
    std::vector<double> theta_values;
    for (const auto& data : surface_data) {
        if (std::find(theta_values.begin(), theta_values.end(), data.theta) == 
            theta_values.end()) {
            theta_values.push_back(data.theta);
        }
    }
    n_theta = theta_values.size();
    
    // Count points per θ (assume constant)
    size_t points_per_theta = 0;
    for (const auto& data : surface_data) {
        if (std::abs(data.theta - theta_values[0]) < 1e-10) {
            points_per_theta++;
        }
    }
    n_phi = points_per_theta;
    
    // Generate faces for θ-φ grid
    if (n_theta > 1 && n_phi > 1) {
        // Map points to grid indices
        std::vector<std::vector<size_t>> grid_indices(n_theta, 
                                                     std::vector<size_t>(n_phi, 0));
        
        // This is simplified - would need proper mapping in production
        
        // Generate triangular faces
        for (size_t i = 0; i < n_theta - 1; i++) {
            for (size_t j = 0; j < n_phi - 1; j++) {
                // Two triangles per quad
                faces.push_back({i*n_phi + j, (i+1)*n_phi + j, i*n_phi + j+1});
                faces.push_back({(i+1)*n_phi + j, (i+1)*n_phi + j+1, i*n_phi + j+1});
            }
        }
    }
    
    return faces;
}

void SurfaceIntegral::interpolate_normals(
    std::vector<MetricData>& surface_data) {
    
    // Compute normals for points that don't have them
    auto computed_normals = compute_surface_normals(surface_data);
    
    for (size_t i = 0; i < surface_data.size(); i++) {
        if (surface_data[i].normal.norm() < 1e-10) {
            surface_data[i].normal = computed_normals[i];
        }
    }
}

} // namespace grid
} // namespace qlt
} // namespace grid 
} // namespace qlt 
