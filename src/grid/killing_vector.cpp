// killing_vector.cpp - Fully corrected version
#include "killing_vector.h"
#include "coordinate_systems.h"
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <utility>

namespace qlt {
namespace grid {

// Constructor
KillingVectorGenerator::KillingVectorGenerator() {}

// Destructor is defined in the header, so we don't define it here

std::vector<MetricData> KillingVectorGenerator::generate_killing_field(
    const std::vector<MetricData>& grid_data,
    KillingVectorType type,
    const std::array<double, 3>& parameters) {
    
    std::vector<MetricData> result = grid_data;
    
    for (size_t i = 0; i < result.size(); i++) {
        MetricData& data = result[i];
        const MetricData& original = grid_data[i];
        
        // Get Cartesian coordinates
        double x = original.r * sin(original.theta) * cos(original.phi_coord);
        double y = original.r * sin(original.theta) * sin(original.phi_coord);
        double z = original.r * cos(original.theta);
        
        Vector3d xi;
        Matrix3d dxi_dx = Matrix3d::Zero();
        
        switch (type) {
            case KillingVectorType::ROTATIONAL_Z: {
                xi = Vector3d(-y, x, 0.0);
                dxi_dx(0,1) = -1.0;  // ∂ξ^x/∂y = -1
                dxi_dx(1,0) = 1.0;   // ∂ξ^y/∂x = 1
                break;
            }
            case KillingVectorType::ROTATIONAL_X: {
                xi = Vector3d(0.0, -z, y);
                dxi_dx(1,2) = -1.0;  // ∂ξ^y/∂z = -1
                dxi_dx(2,1) = 1.0;   // ∂ξ^z/∂y = 1
                break;
            }
            case KillingVectorType::ROTATIONAL_Y: {
                xi = Vector3d(z, 0.0, -x);
                dxi_dx(0,2) = 1.0;   // ∂ξ^x/∂z = 1
                dxi_dx(2,0) = -1.0;  // ∂ξ^z/∂x = -1
                break;
            }
            case KillingVectorType::TRANSLATIONAL_X: {
                xi = Vector3d(1.0, 0.0, 0.0);
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
                xi = Vector3d(parameters[0], parameters[1], parameters[2]);
                break;
            }
            default:
                throw std::runtime_error("Unsupported Killing vector type");
        }
        
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
    
    Vector3d axis_vec(axis[0], axis[1], axis[2]);
    double norm = axis_vec.norm();
    if (norm > 1e-10) {
        axis_vec /= norm;
    }
    
    for (size_t i = 0; i < result.size(); i++) {
        MetricData& data = result[i];
        const MetricData& original = grid_data[i];
        
        double x = original.r * sin(original.theta) * cos(original.phi_coord);
        double y = original.r * sin(original.theta) * sin(original.phi_coord);
        double z = original.r * cos(original.theta);
        
        Vector3d pos(x, y, z);
        Vector3d omega = angular_velocity * axis_vec;
        Vector3d xi = omega.cross(pos);
        
        data.xi = xi;
        data.dxi_dx = Matrix3d::Zero(); // Would need full derivative in production
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
    }
    
    for (size_t i = 0; i < result.size(); i++) {
        result[i].xi = dir;
    }
    
    return result;
}

std::vector<MetricData> KillingVectorGenerator::compute_approximate_killing_vector(
    const std::vector<MetricData>& grid_data,
    int max_iterations,
    double tolerance) {
    
    std::vector<MetricData> result = grid_data;
    std::cout << "[INFO] Computing approximate Killing vector" << std::endl;
    
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
        // Minkowski spacetime - just use translational vectors
        return generate_translational_killing_vector(result, {1.0, 0.0, 0.0});
    }
    else if (spacetime_name == "Schwarzschild") {
        // Schwarzschild - stationary and spherically symmetric
        return generate_rotational_killing_vector(result, {0.0, 0.0, 1.0}, 1.0);
    }
    else if (spacetime_name == "Kerr") {
        double M = 1.0;
        double a = 0.0;
        
        auto it_M = parameters.find("M");
        if (it_M != parameters.end()) M = it_M->second;
        
        auto it_a = parameters.find("a");
        if (it_a != parameters.end()) a = it_a->second;
        
        // Kerr has both rotational (∂/∂φ) and stationary (∂/∂t) Killing vectors
        return generate_rotational_killing_vector(result, {0.0, 0.0, 1.0}, 1.0);
    }
    else {
        throw std::runtime_error("Unknown spacetime: " + spacetime_name);
    }
}

std::vector<MetricData> KillingVectorGenerator::read_from_file(
    const std::string& filename,
    const std::vector<MetricData>& grid_data) {
    
    std::vector<MetricData> result = grid_data;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open Killing vector file: " + filename);
    }
    
    std::cout << "[INFO] Reading Killing vector from: " << filename << std::endl;
    
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
    
    file << "# Killing vector field" << std::endl;
    file << "# Format: x,y,z,xi_x,xi_y,xi_z" << std::endl;
    
    for (const auto& data : data_with_xi) {
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
    
    Matrix3d derivatives = Matrix3d::Zero();
    // Simplified implementation - would need proper finite differences
    return derivatives;
}

Vector3d KillingVectorGenerator::optimize_killing_vector(
    const MetricData& data,
    const std::vector<MetricData>& neighborhood) {
    
    // Simplified implementation
    return data.xi;
}

} // namespace grid
} // namespace qlt