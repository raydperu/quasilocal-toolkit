// killing_vector.cpp 
#include "killing_vector.h" 
 
namespace qlt { 
namespace grid { 
 
#include "killing_vector.h"
#include "coordinate_systems.h"
#include <algorithm>
#include <cmath>
#include <fstream>

namespace qlt {
namespace grid {

KillingVectorGenerator::KillingVectorGenerator() = default;

std::vector<MetricData> KillingVectorGenerator::generate_killing_field(
    const std::vector<MetricData>& grid_data,
    KillingVectorType type,
    const std::array<double, 3>& parameters) {
    
    std::vector<MetricData> result = grid_data;
    
    for (size_t i = 0; i < result.size(); i++) {
        MetricData& data = result[i];
        const auto& original = grid_data[i];
        
        Vector3d xi;
        Matrix3d dxi_dx = Matrix3d::Zero();
        
        switch (type) {
            case KillingVectorType::ROTATIONAL_Z: {
                // ξ = ∂/∂φ = (-y, x, 0) in Cartesian coordinates
                // Convert from spherical to Cartesian if needed
                double x = original.r * sin(original.theta) * cos(original.phi_coord);
                double y = original.r * sin(original.theta) * sin(original.phi_coord);
                
                xi = Vector3d(-y, x, 0.0);
                
                // Derivatives: ∂ξ^x/∂x = 0, ∂ξ^x/∂y = -1, ∂ξ^x/∂z = 0
                //              ∂ξ^y/∂x = 1, ∂ξ^y/∂y = 0, ∂ξ^y/∂z = 0
                dxi_dx(0, 1) = -1.0;  // ∂ξ^x/∂y
                dxi_dx(1, 0) = 1.0;   // ∂ξ^y/∂x
                break;
            }
            
            case KillingVectorType::ROTATIONAL_X: {
                // Rotation about x-axis: ξ = (0, -z, y)
                double y = original.r * sin(original.theta) * sin(original.phi_coord);
                double z = original.r * cos(original.theta);
                
                xi = Vector3d(0.0, -z, y);
                dxi_dx(1, 2) = -1.0;  // ∂ξ^y/∂z
                dxi_dx(2, 1) = 1.0;   // ∂ξ^z/∂y
                break;
            }
            
            case KillingVectorType::ROTATIONAL_Y: {
                // Rotation about y-axis: ξ = (z, 0, -x)
                double x = original.r * sin(original.theta) * cos(original.phi_coord);
                double z = original.r * cos(original.theta);
                
                xi = Vector3d(z, 0.0, -x);
                dxi_dx(0, 2) = 1.0;   // ∂ξ^x/∂z
                dxi_dx(2, 0) = -1.0;  // ∂ξ^z/∂x
                break;
            }
            
            case KillingVectorType::TRANSLATIONAL_X: {
                // ξ = ∂/∂x = (1, 0, 0)
                xi = Vector3d(1.0, 0.0, 0.0);
                // Derivatives are all zero
                break;
            }
            
            case KillingVectorType::TRANSLATIONAL_Y: {
                xi = Vector3d(0.0, 1.0, 0.0);
                break;
            }
            
            case KillingVectorType::TRANSLATIONAL_Z: {
                xi = Vector3d(0.0, 0.0, 1.0);
                break;
            }
            
            case KillingVectorType::CUSTOM: {
                // User-provided parameters as (vx, vy, vz)
                xi = Vector3d(parameters[0], parameters[1], parameters[2]);
                break;
            }
            
            default:
                throw std::runtime_error("Unsupported Killing vector type");
        }
        
        // Store in data structure
        data.xi = xi;
        data.dxi_dx = dxi_dx;
    }
    
    return result;
}

std::vector<MetricData> KillingVectorGenerator::generate_rotational_killing_vector(
    const std::vector<MetricData>& grid_data,
    const std::array<double, 3>& axis,
    double angular_velocity) {
    
    std::vector<MetricData> result = grid_data;
    
    // Normalize axis
    Vector3d axis_vec(axis[0], axis[1], axis[2]);
    double norm = axis_vec.norm();
    if (norm > 1e-10) {
        axis_vec /= norm;
    } else {
        axis_vec = Vector3d(0.0, 0.0, 1.0);  // Default to z-axis
    }
    
    for (size_t i = 0; i < result.size(); i++) {
        const auto& original = grid_data[i];
        
        // Position vector in Cartesian coordinates
        double x = original.r * sin(original.theta) * cos(original.phi_coord);
        double y = original.r * sin(original.theta) * sin(original.phi_coord);
        double z = original.r * cos(original.theta);
        
        Vector3d pos(x, y, z);
        
        // Rotational Killing vector: ξ = Ω × r
        // where Ω = angular_velocity * axis
        Vector3d omega = angular_velocity * axis_vec;
        Vector3d xi = omega.cross(pos);
        
        // Compute derivatives: ∇ξ is antisymmetric for rotation
        Matrix3d dxi_dx = Matrix3d::Zero();
        
        // For Ω = (0,0,ω): dξ^x/dy = -ω, dξ^y/dx = ω
        // General case: dξ_i/dx_j = ε_{ijk} Ω^k
        dxi_dx(0, 1) = -omega[2];  // ∂ξ^x/∂y
        dxi_dx(0, 2) = omega[1];   // ∂ξ^x/∂z
        dxi_dx(1, 0) = omega[2];   // ∂ξ^y/∂x
        dxi_dx(1, 2) = -omega[0];  // ∂ξ^y/∂z
        dxi_dx(2, 0) = -omega[1];  // ∂ξ^z/∂x
        dxi_dx(2, 1) = omega[0];   // ∂ξ^z/∂y
        
        result[i].xi = xi;
        result[i].dxi_dx = dxi_dx;
    }
    
    return result;
}

std::vector<MetricData> KillingVectorGenerator::generate_translational_killing_vector(
    const std::vector<MetricData>& grid_data,
    const std::array<double, 3>& direction) {
    
    std::vector<MetricData> result = grid_data;
    
    Vector3d dir(direction[0], direction[1], direction[2]);
    double norm = dir.norm();
    if (norm > 1e-10) {
        dir /= norm;
    } else {
        dir = Vector3d(1.0, 0.0, 0.0);  // Default to x-direction
    }
    
    for (size_t i = 0; i < result.size(); i++) {
        result[i].xi = dir;
        result[i].dxi_dx = Matrix3d::Zero();  // Constant vector field
    }
    
    return result;
}

std::vector<MetricData> KillingVectorGenerator::compute_approximate_killing_vector(
    const std::vector<MetricData>& grid_data,
    int max_iterations,
    double tolerance) {
    
    std::cout << "[INFO] Computing approximate Killing vector" << std::endl;
    
    // This is a simplified implementation
    // In production, would use eigen-decomposition of Killing operator
    
    std::vector<MetricData> result = grid_data;
    
    // Initialize with rotational guess
    result = generate_rotational_killing_vector(grid_data, {0.0, 0.0, 1.0});
    
    // Simple optimization loop (conceptual)
    for (int iter = 0; iter < max_iterations; iter++) {
        double max_violation = 0.0;
        double rms_violation = 0.0;
        
        check_killing_equation(result, max_violation, rms_violation);
        
        std::cout << "[INFO] Iteration " << iter 
                  << ": max violation = " << max_violation
                  << ", RMS = " << rms_violation << std::endl;
        
        if (rms_violation < tolerance) {
            std::cout << "[INFO] Converged after " << iter << " iterations" << std::endl;
            break;
        }
        
        // In full implementation, would adjust ξ to reduce violation
        // This would involve solving a linear system
    }
    
    return result;
}

std::pair<double, double> KillingVectorGenerator::check_killing_equation(
    const std::vector<MetricData>& data_with_xi,
    double& max_violation,
    double& rms_violation) {
    
    max_violation = 0.0;
    rms_violation = 0.0;
    int count = 0;
    
    for (const auto& data : data_with_xi) {
        // Compute Killing equation violation: ∇_i ξ_j + ∇_j ξ_i
        // For exact Killing vector, this should be zero
        
        // Simplified check using provided derivatives
        // In production, would compute covariant derivative properly
        
        Matrix3d killing_tensor = data.dxi_dx + data.dxi_dx.transpose();
        double norm = killing_tensor.norm();
        
        max_violation = std::max(max_violation, norm);
        rms_violation += norm * norm;
        count++;
    }
    
    if (count > 0) {
        rms_violation = std::sqrt(rms_violation / count);
    }
    
    return std::make_pair(max_violation, rms_violation);
}

std::vector<Matrix3d> KillingVectorGenerator::compute_lie_derivative(
    const std::vector<MetricData>& data) {
    
    std::vector<Matrix3d> lie_derivs;
    lie_derivs.reserve(data.size());
    
    for (const auto& d : data) {
        // Lie derivative of metric along ξ: £_ξ g_{ij} = ∇_i ξ_j + ∇_j ξ_i
        // For Killing vector, this should be zero
        
        Matrix3d lie_deriv = d.dxi_dx + d.dxi_dx.transpose();
        lie_derivs.push_back(lie_deriv);
    }
    
    return lie_derivs;
}

std::vector<MetricData> KillingVectorGenerator::generate_for_spacetime(
    const std::vector<MetricData>& grid_data,
    const std::string& spacetime_name,
    const std::map<std::string, double>& parameters) {
    
    std::vector<MetricData> result = grid_data;
    
    if (spacetime_name == "Minkowski" || spacetime_name == "Flat") {
        // Minkowski has 10 Killing vectors (Poincaré group)
        // Return rotational about z by default
        return generate_rotational_killing_vector(grid_data);
    }
    else if (spacetime_name == "Schwarzschild") {
        // Schwarzschild is static and spherically symmetric
        // Has timelike Killing vector ∂/∂t and 3 rotational Killing vectors
        return generate_rotational_killing_vector(grid_data);
    }
    else if (spacetime_name == "Kerr") {
        double M = 1.0;
        double a = 0.0;
        
        auto it_M = parameters.find("M");
        if (it_M != parameters.end()) M = it_M->second;
        
        auto it_a = parameters.find("a");
        if (it_a != parameters.end()) a = it_a->second;
        
        // Kerr has 2 Killing vectors: timelike ∂/∂t and rotational ∂/∂φ
        // The rotational one is what we use for angular momentum
        return generate_rotational_killing_vector(grid_data);
    }
    else {
        throw std::runtime_error("Unknown spacetime: " + spacetime_name);
    }
}

std::vector<MetricData> KillingVectorGenerator::read_from_file(
    const std::string& filename,
    const std::vector<MetricData>& grid_data) {
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open Killing vector file: " + filename);
    }
    
    std::cout << "[INFO] Reading Killing vector from: " << filename << std::endl;
    
    std::vector<MetricData> result = grid_data;
    
    // Simple CSV format: x,y,z,xi_x,xi_y,xi_z
    std::string line;
    size_t line_num = 0;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        if (line_num >= result.size()) {
            std::cerr << "[WARNING] More lines in file than grid points" << std::endl;
            break;
        }
        
        std::istringstream iss(line);
        std::string token;
        std::vector<double> values;
        
        while (std::getline(iss, token, ',')) {
            try {
                values.push_back(std::stod(token));
            } catch (const std::exception& e) {
                std::cerr << "[WARNING] Could not parse value on line " 
                          << line_num << ": " << token << std::endl;
                values.push_back(0.0);
            }
        }
        
        if (values.size() >= 6) {
            // Assume format: x,y,z,xi_x,xi_y,xi_z
            result[line_num].xi = Vector3d(values[3], values[4], values[5]);
        }
        
        line_num++;
    }
    
    std::cout << "[INFO] Read " << line_num << " Killing vectors" << std::endl;
    
    return result;
}

void KillingVectorGenerator::save_to_file(
    const std::vector<MetricData>& data_with_xi,
    const std::string& filename) {
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    // Write header
    file << "# Killing vector field" << std::endl;
    file << "# Format: x,y,z,xi_x,xi_y,xi_z" << std::endl;
    
    for (const auto& data : data_with_xi) {
        // Convert to Cartesian coordinates
        double x = data.r * sin(data.theta) * cos(data.phi_coord);
        double y = data.r * sin(data.theta) * sin(data.phi_coord);
        double z = data.r * cos(data.theta);
        
        file << std::scientific << std::setprecision(15)
             << x << "," << y << "," << z << ","
             << data.xi[0] << "," << data.xi[1] << "," << data.xi[2]
             << std::endl;
    }
    
    std::cout << "[INFO] Saved " << data_with_xi.size() 
              << " Killing vectors to: " << filename << std::endl;
}

Matrix3d KillingVectorGenerator::compute_killing_derivatives(
    const MetricData& data,
    const std::vector<MetricData>& neighborhood) {
    
    // Finite difference computation of ∇ξ
    // In production, would use proper stencil
    
    return data.dxi_dx;  // Use provided for now
}

Vector3d KillingVectorGenerator::optimize_killing_vector(
    const MetricData& data,
    const std::vector<MetricData>& neighborhood) {
    
    // Optimization to find best Killing vector at point
    // In production, would solve minimization problem
    
    return data.xi;  // Return current for now
}

} // namespace grid
} // namespace qlt
} // namespace grid 
} // namespace qlt 
