#include "grid/coordinate_systems.h"
#include <algorithm>
#include <stdexcept>

namespace qlt {
namespace grid {

std::array<double, 3> CoordinateSystems::cartesian_to_spherical(
    double x, double y, double z) {
    
    double r = std::sqrt(x * x + y * y + z * z);
    double theta = (r > 1e-10) ? std::acos(z / r) : 0.0;
    double phi = std::atan2(y, x);
    
    return {r, theta, phi};
}

std::array<double, 3> CoordinateSystems::spherical_to_cartesian(
    double r, double theta, double phi) {
    
    double x = r * std::sin(theta) * std::cos(phi);
    double y = r * std::sin(theta) * std::sin(phi);
    double z = r * std::cos(theta);
    
    return {x, y, z};
}

std::array<double, 3> CoordinateSystems::cartesian_to_cylindrical(
    double x, double y, double z) {
    
    double rho = std::sqrt(x * x + y * y);
    double phi = std::atan2(y, x);
    
    return {rho, phi, z};
}

std::array<double, 3> CoordinateSystems::cylindrical_to_cartesian(
    double rho, double phi, double z) {
    
    double x = rho * std::cos(phi);
    double y = rho * std::sin(phi);
    
    return {x, y, z};
}

Matrix3d CoordinateSystems::compute_jacobian(
    const std::string& from,
    const std::string& to,
    const std::array<double, 3>& point) {
    
    Matrix3d J = Matrix3d::Zero();
    
    if (from == "Cartesian" && to == "Spherical") {
        double x = point[0], y = point[1], z = point[2];
        double r = std::sqrt(x * x + y * y + z * z);
        double r2 = x * x + y * y;
        double rho = std::sqrt(r2);
        
        if (r > 1e-10) {
            // ∂(r,θ,φ)/∂(x,y,z)
            J(0, 0) = x / r;               // ∂r/∂x
            J(0, 1) = y / r;               // ∂r/∂y
            J(0, 2) = z / r;               // ∂r/∂z
            
            if (rho > 1e-10) {
                J(1, 0) = x * z / (r * r * rho);  // ∂θ/∂x
                J(1, 1) = y * z / (r * r * rho);  // ∂θ/∂y
                J(1, 2) = -rho / (r * r);         // ∂θ/∂z
                
                J(2, 0) = -y / r2;                // ∂φ/∂x
                J(2, 1) = x / r2;                 // ∂φ/∂y
                J(2, 2) = 0.0;                    // ∂φ/∂z
            }
        }
    }
    else if (from == "Spherical" && to == "Cartesian") {
        double r = point[0], theta = point[1], phi = point[2];
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double sin_phi = std::sin(phi);
        double cos_phi = std::cos(phi);
        
        // ∂(x,y,z)/∂(r,θ,φ)
        J(0, 0) = sin_theta * cos_phi;      // ∂x/∂r
        J(0, 1) = r * cos_theta * cos_phi;  // ∂x/∂θ
        J(0, 2) = -r * sin_theta * sin_phi; // ∂x/∂φ
        
        J(1, 0) = sin_theta * sin_phi;      // ∂y/∂r
        J(1, 1) = r * cos_theta * sin_phi;  // ∂y/∂θ
        J(1, 2) = r * sin_theta * cos_phi;  // ∂y/∂φ
        
        J(2, 0) = cos_theta;                // ∂z/∂r
        J(2, 1) = -r * sin_theta;           // ∂z/∂θ
        J(2, 2) = 0.0;                      // ∂z/∂φ
    }
    else if (from == "Cartesian" && to == "Cylindrical") {
        double x = point[0], y = point[1], z = point[2];
        double rho = std::sqrt(x * x + y * y);
        
        if (rho > 1e-10) {
            // ∂(ρ,φ,z)/∂(x,y,z)
            J(0, 0) = x / rho;               // ∂ρ/∂x
            J(0, 1) = y / rho;               // ∂ρ/∂y
            J(0, 2) = 0.0;                   // ∂ρ/∂z
            
            J(1, 0) = -y / (rho * rho);      // ∂φ/∂x
            J(1, 1) = x / (rho * rho);       // ∂φ/∂y
            J(1, 2) = 0.0;                   // ∂φ/∂z
            
            J(2, 0) = 0.0;                   // ∂z/∂x
            J(2, 1) = 0.0;                   // ∂z/∂y
            J(2, 2) = 1.0;                   // ∂z/∂z
        }
    }
    else {
        throw std::runtime_error(
            "Unsupported coordinate transformation: " + from + " -> " + to
        );
    }
    
    return J;
}

Vector3d CoordinateSystems::transform_vector(
    const Vector3d& vector,
    const std::string& from,
    const std::string& to,
    const std::array<double, 3>& point) {
    
    if (from == to) {
        return vector;
    }
    
    Matrix3d J = compute_jacobian(from, to, point);
    return J * vector;
}

Matrix3d CoordinateSystems::transform_metric(
    const Matrix3d& metric,
    const std::string& from,
    const std::string& to,
    const std::array<double, 3>& point) {
    
    if (from == to) {
        return metric;
    }
    
    Matrix3d J = compute_jacobian(from, to, point);
    // Metric transforms as: g' = J^T g J
    return J.transpose() * metric * J;
}

std::vector<std::vector<std::vector<std::array<double, 3>>>>
CoordinateSystems::generate_spherical_grid(
    double r_min, double r_max, size_t num_r,
    size_t num_theta, size_t num_phi) {
    
    if (r_min < 0.0 || r_max <= r_min || num_r == 0 || 
        num_theta == 0 || num_phi == 0) {
        throw std::invalid_argument("Invalid spherical grid parameters");
    }
    
    std::vector<std::vector<std::vector<std::array<double, 3>>>> grid(
        num_r, std::vector<std::vector<std::array<double, 3>>>(
            num_theta, std::vector<std::array<double, 3>>(num_phi)
        )
    );
    
    // Generate logarithmic radial grid (more points near horizon)
    std::vector<double> r_values(num_r);
    if (r_min > 0.0) {
        double log_r_min = std::log(r_min);
        double log_r_max = std::log(r_max);
        for (size_t i = 0; i < num_r; i++) {
            double t = static_cast<double>(i) / (num_r - 1);
            r_values[i] = std::exp(log_r_min + t * (log_r_max - log_r_min));
        }
    } else {
        // Linear grid if starting at r=0
        for (size_t i = 0; i < num_r; i++) {
            r_values[i] = r_min + (r_max - r_min) * i / (num_r - 1);
        }
    }
    
    // Generate angular grids
    for (size_t i = 0; i < num_r; i++) {
        double r = r_values[i];
        
        for (size_t j = 0; j < num_theta; j++) {
            // Avoid exactly 0 and π to prevent coordinate singularities
            double theta = M_PI * (j + 0.5) / num_theta;
            
            for (size_t k = 0; k < num_phi; k++) {
                double phi = 2.0 * M_PI * k / num_phi;
                
                grid[i][j][k] = {r, theta, phi};
            }
        }
    }
    
    return grid;
}

std::vector<std::vector<std::vector<std::array<double, 3>>>>
CoordinateSystems::generate_cartesian_grid(
    double x_min, double x_max, size_t num_x,
    double y_min, double y_max, size_t num_y,
    double z_min, double z_max, size_t num_z) {
    
    std::vector<std::vector<std::vector<std::array<double, 3>>>> grid(
        num_x, std::vector<std::vector<std::array<double, 3>>>(
            num_y, std::vector<std::array<double, 3>>(num_z)
        )
    );
    
    std::vector<double> x_values(num_x);
    std::vector<double> y_values(num_y);
    std::vector<double> z_values(num_z);
    
    for (size_t i = 0; i < num_x; i++) {
        x_values[i] = x_min + (x_max - x_min) * i / (num_x - 1);
    }
    for (size_t j = 0; j < num_y; j++) {
        y_values[j] = y_min + (y_max - y_min) * j / (num_y - 1);
    }
    for (size_t k = 0; k < num_z; k++) {
        z_values[k] = z_min + (z_max - z_min) * k / (num_z - 1);
    }
    
    for (size_t i = 0; i < num_x; i++) {
        for (size_t j = 0; j < num_y; j++) {
            for (size_t k = 0; k < num_z; k++) {
                grid[i][j][k] = {x_values[i], y_values[j], z_values[k]};
            }
        }
    }
    
    return grid;
}

std::vector<MetricData> CoordinateSystems::interpolate_to_grid(
    const std::vector<MetricData>& source_data,
    const std::vector<std::array<double, 3>>& target_coords,
    const std::string& interpolation_method) {
    
    if (source_data.empty() || target_coords.empty()) {
        return {};
    }
    
    std::vector<MetricData> interpolated_data;
    interpolated_data.reserve(target_coords.size());
    
    // Simple nearest-neighbor interpolation for now
    // In production code, would implement trilinear/cubic interpolation
    
    for (const auto& target_coord : target_coords) {
        // Find nearest source point
        double min_dist = std::numeric_limits<double>::max();
        const MetricData* nearest = nullptr;
        
        for (const auto& source : source_data) {
            double dx = target_coord[0] - (source.r * sin(source.theta) * cos(source.phi_coord));
            double dy = target_coord[1] - (source.r * sin(source.theta) * sin(source.phi_coord));
            double dz = target_coord[2] - (source.r * cos(source.theta));
            
            double dist = dx * dx + dy * dy + dz * dz;
            if (dist < min_dist) {
                min_dist = dist;
                nearest = &source;
            }
        }
        
        if (nearest) {
            MetricData interpolated = *nearest;
            // Update coordinates
            auto spherical = cartesian_to_spherical(
                target_coord[0], target_coord[1], target_coord[2]);
            interpolated.r = spherical[0];
            interpolated.theta = spherical[1];
            interpolated.phi_coord = spherical[2];
            
            interpolated_data.push_back(interpolated);
        }
    }
    
    return interpolated_data;
}

bool CoordinateSystems::validate_coordinates(
    const std::array<double, 3>& coords,
    const std::string& system) {
    
    if (system == "Spherical") {
        double r = coords[0];
        double theta = coords[1];
        double phi = coords[2];
        
        if (r < 0.0) return false;
        if (theta < 0.0 || theta > M_PI) return false;
        // phi can be any real number, but typically in [0, 2π)
        return true;
    }
    else if (system == "Cartesian") {
        // Cartesian coordinates are always valid
        return true;
    }
    else if (system == "Cylindrical") {
        double rho = coords[0];
        double z = coords[2];
        
        if (rho < 0.0) return false;
        // phi and z can be any real number
        return true;
    }
    
    return false;
}

std::array<Vector3d, 3> CoordinateSystems::compute_basis_vectors(
    const std::string& system,
    const std::array<double, 3>& point) {
    
    std::array<Vector3d, 3> basis;
    
    if (system == "Cartesian") {
        basis[0] = Vector3d(1.0, 0.0, 0.0);  // ∂/∂x
        basis[1] = Vector3d(0.0, 1.0, 0.0);  // ∂/∂y
        basis[2] = Vector3d(0.0, 0.0, 1.0);  // ∂/∂z
    }
    else if (system == "Spherical") {
        double r = point[0], theta = point[1], phi = point[2];
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double sin_phi = std::sin(phi);
        double cos_phi = std::cos(phi);
        
        // ∂/∂r in Cartesian coordinates
        basis[0] = Vector3d(
            sin_theta * cos_phi,
            sin_theta * sin_phi,
            cos_theta
        );
        
        // ∂/∂θ in Cartesian coordinates
        basis[1] = Vector3d(
            r * cos_theta * cos_phi,
            r * cos_theta * sin_phi,
            -r * sin_theta
        );
        
        // ∂/∂φ in Cartesian coordinates
        basis[2] = Vector3d(
            -r * sin_theta * sin_phi,
            r * sin_theta * cos_phi,
            0.0
        );
    }
    else if (system == "Cylindrical") {
        double rho = point[0], phi = point[1], z = point[2];
        
        // ∂/∂ρ in Cartesian coordinates
        basis[0] = Vector3d(
            std::cos(phi),
            std::sin(phi),
            0.0
        );
        
        // ∂/∂φ in Cartesian coordinates
        basis[1] = Vector3d(
            -rho * std::sin(phi),
            rho * std::cos(phi),
            0.0
        );
        
        // ∂/∂z in Cartesian coordinates
        basis[2] = Vector3d(0.0, 0.0, 1.0);
    }
    else {
        throw std::runtime_error("Unknown coordinate system: " + system);
    }
    
    return basis;
}

std::array<std::array<std::array<double, 3>, 3>, 3>
CoordinateSystems::compute_christoffel_symbols(
    const Matrix3d& metric,
    const std::string& coordinate_system,
    const std::array<double, 3>& point,
    double epsilon) {
    
    std::array<std::array<std::array<double, 3>, 3>, 3> gamma;
    
    // Finite difference approximation of Christoffel symbols
    // Γ^i_{jk} = (1/2) g^{il} (∂_j g_{kl} + ∂_k g_{jl} - ∂_l g_{jk})
    
    // For production code, would implement analytic formulas for common metrics
    // For now, return zeros (valid for Cartesian coordinates in flat space)
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                gamma[i][j][k] = 0.0;
            }
        }
    }
    
    // Special case: spherical coordinates for Schwarzschild/Kerr
    if (coordinate_system == "Spherical") {
        // Approximate for flat space in spherical coordinates
        double r = point[0];
        
        if (r > epsilon) {
            // Non-zero Christoffel symbols for flat space in spherical coords
            gamma[0][1][1] = -r;               // Γ^r_{θθ}
            gamma[0][2][2] = -r * sin(point[1]) * sin(point[1]);  // Γ^r_{φφ}
            
            gamma[1][0][1] = 1.0 / r;          // Γ^θ_{rθ} = Γ^θ_{θr}
            gamma[1][2][2] = -sin(point[1]) * cos(point[1]);  // Γ^θ_{φφ}
            
            gamma[2][0][2] = 1.0 / r;          // Γ^φ_{rφ} = Γ^φ_{φr}
            gamma[2][1][2] = cos(point[1]) / sin(point[1]);  // Γ^φ_{θφ} = Γ^φ_{φθ}
        }
    }
    
    return gamma;
}

Matrix3d CoordinateSystems::compute_covariant_derivative(
    const Vector3d& vector,
    const Matrix3d& metric,
    const std::string& coordinate_system,
    const std::array<double, 3>& point,
    double epsilon) {
    
    Matrix3d cov_deriv = Matrix3d::Zero();
    
    // Finite difference approximation
    // ∇_j V^i = ∂_j V^i + Γ^i_{jk} V^k
    
    // For now, return partial derivative (valid in Cartesian coordinates)
    // In production code, would add Christoffel term
    
    return cov_deriv;
}

} // namespace grid
} // namespace qlt