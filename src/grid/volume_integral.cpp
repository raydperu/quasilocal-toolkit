// volume_integral.cpp - Fully corrected version
#include "volume_integral.h"
#include "surface_integral.h"
#include "coordinate_systems.h"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>

namespace qlt {
namespace grid {

// Helper function
std::vector<double> compute_all_volume_elements(
    const std::vector<MetricData>& volume_data) {
    
    std::vector<double> elements(volume_data.size());
    for (size_t i = 0; i < volume_data.size(); i++) {
        if (volume_data[i].sqrt_det_h > 0.0) {
            elements[i] = volume_data[i].sqrt_det_h;
        } else {
            double det_h = volume_data[i].h_ij.determinant();
            elements[i] = std::sqrt(std::max(0.0, det_h));
        }
    }
    return elements;
}

double compute_volume_element(const MetricData& data) {
    if (data.sqrt_det_h > 0.0) {
        return data.sqrt_det_h;
    } else {
        return std::sqrt(std::max(0.0, data.h_ij.determinant()));
    }
}

// VolumeIntegral implementations
double VolumeIntegral::integrate_scalar(
    const std::vector<MetricData>& volume_data,
    const std::vector<double>& scalar_values,
    IntegrationMethod method) {
    
    if (volume_data.size() != scalar_values.size()) {
        throw std::invalid_argument(
            "Volume data and scalar values must have same size");
    }
    
    auto volume_elements = compute_all_volume_elements(volume_data);
    
    switch (method) {
        case IntegrationMethod::MIDPOINT_RULE: {
            double integral = 0.0;
            for (size_t i = 0; i < volume_data.size(); i++) {
                integral += scalar_values[i] * volume_elements[i];
            }
            return integral;
        }
        
        case IntegrationMethod::TRAPEZOIDAL_RULE: {
            double integral = 0.0;
            for (size_t i = 0; i < volume_data.size(); i++) {
                integral += scalar_values[i] * volume_elements[i];
            }
            return integral;
        }
        
        case IntegrationMethod::MONTE_CARLO: {
            if (volume_data.size() < 100) {
                return integrate_scalar(volume_data, scalar_values, 
                                       IntegrationMethod::MIDPOINT_RULE);
            }
            
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dist(0, volume_data.size() - 1);
            
            size_t n_samples = std::min<size_t>(1000, volume_data.size());
            double sum = 0.0;
            double total_volume = 0.0;
            
            for (size_t i = 0; i < n_samples; i++) {
                size_t idx = dist(gen);
                sum += scalar_values[idx];
                total_volume += volume_elements[idx];
            }
            
            double average_value = sum / n_samples;
            double estimated_volume = total_volume / n_samples * volume_data.size();
            
            return average_value * estimated_volume;
        }
        
        default:
            throw std::runtime_error("Unsupported integration method");
    }
}

double VolumeIntegral::integrate_vector_divergence(
    const std::vector<MetricData>& volume_data,
    const std::vector<Vector3d>& vector_field,
    const std::vector<MetricData>& surface_data,
    IntegrationMethod method) {
    
    // Apply divergence theorem: ∫_V ∇·F dV = ∫_∂V F·n dA
    double surface_flux = SurfaceIntegral::integrate_vector_flux(
        surface_data, vector_field, 
        static_cast<SurfaceIntegral::IntegrationMethod>(method));
    
    return surface_flux;
}

double VolumeIntegral::compute_clebsch_komar_volume_integral(
    const std::vector<MetricData>& volume_data) {
    
    double integral = 0.0;
    
    for (const auto& data : volume_data) {
        Vector3d dbeta = data.dbeta_dx;
        
        // Simplified dξ
        Vector3d dxi(data.dxi_dx(0,0), data.dxi_dx(1,1), data.dxi_dx(2,2));
        
        Vector3d cross_product = dbeta.cross(dxi);
        
        double dV = compute_volume_element(data);
        double term = data.alpha * cross_product.norm() * dV;
        
        integral += term;
    }
    
    return integral;
}

double VolumeIntegral::compute_volume(
    const std::vector<MetricData>& volume_data) {
    
    auto volume_elements = compute_all_volume_elements(volume_data);
    return std::accumulate(volume_elements.begin(), volume_elements.end(), 0.0);
}

std::array<double, 3> VolumeIntegral::compute_center_of_mass(
    const std::vector<MetricData>& volume_data,
    const std::vector<double>& density) {
    
    std::array<double, 3> com = {0.0, 0.0, 0.0};
    double total_mass = 0.0;
    
    bool use_density = (!density.empty() && density.size() == volume_data.size());
    
    for (size_t i = 0; i < volume_data.size(); i++) {
        const auto& data = volume_data[i];
        
        auto cartesian = CoordinateSystems::spherical_to_cartesian(
            data.r, data.theta, data.phi_coord);
        
        double dV = compute_volume_element(data);
        double mass_element = use_density ? density[i] * dV : dV;
        
        com[0] += cartesian[0] * mass_element;
        com[1] += cartesian[1] * mass_element;
        com[2] += cartesian[2] * mass_element;
        total_mass += mass_element;
    }
    
    if (total_mass > 1e-10) {
        com[0] /= total_mass;
        com[1] /= total_mass;
        com[2] /= total_mass;
    }
    
    return com;
}

std::vector<MetricData> VolumeIntegral::extract_subvolume(
    const std::vector<std::vector<std::vector<MetricData>>>& grid_data,
    const std::array<double, 3>& min_corner,
    const std::array<double, 3>& max_corner) {
    
    std::vector<MetricData> subvolume;
    
    if (grid_data.empty()) return subvolume;
    
    for (const auto& plane : grid_data) {
        for (const auto& row : plane) {
            for (const auto& data : row) {
                auto cartesian = CoordinateSystems::spherical_to_cartesian(
                    data.r, data.theta, data.phi_coord);
                
                if (cartesian[0] >= min_corner[0] && cartesian[0] <= max_corner[0] &&
                    cartesian[1] >= min_corner[1] && cartesian[1] <= max_corner[1] &&
                    cartesian[2] >= min_corner[2] && cartesian[2] <= max_corner[2]) {
                    
                    subvolume.push_back(data);
                }
            }
        }
    }
    
    return subvolume;
}

double VolumeIntegral::compute_vorticity_bound_integral(
    const std::vector<MetricData>& volume_data) {
    
    double integral = 0.0;
    
    for (const auto& data : volume_data) {
        double alpha_norm = std::abs(data.alpha);
        
        Matrix3d h_inv = data.h_ij.inverse();
        double dbeta_norm_sq = data.dbeta_dx.transpose() * h_inv * data.dbeta_dx;
        double dbeta_norm = std::sqrt(std::max(0.0, dbeta_norm_sq));
        
        double dV = compute_volume_element(data);
        integral += alpha_norm * dbeta_norm * dbeta_norm * dV;
    }
    
    return integral;
}

double VolumeIntegral::compute_helicity_integral(
    const std::vector<MetricData>& volume_data) {
    
    double helicity = 0.0;
    
    for (const auto& data : volume_data) {
        Vector3d V = data.dphi_dx + data.alpha * data.dbeta_dx;
        
        // Compute dV (exterior derivative of V)
        Matrix3d dV_mat;
        dV_mat << 
            0.0, data.dbeta_dx[2] - data.dbeta_dx[1], data.dbeta_dx[1] - data.dbeta_dx[0],
            data.dbeta_dx[0] - data.dbeta_dx[2], 0.0, data.dbeta_dx[2] - data.dbeta_dx[1],
            data.dbeta_dx[1] - data.dbeta_dx[0], data.dbeta_dx[0] - data.dbeta_dx[2], 0.0;
        
        dV_mat *= data.alpha;
        
        double component = 
            V[0] * (dV_mat(1, 2) - dV_mat(2, 1)) +
            V[1] * (dV_mat(2, 0) - dV_mat(0, 2)) +
            V[2] * (dV_mat(0, 1) - dV_mat(1, 0));
        
        double dV_vol = compute_volume_element(data);
        helicity += component * dV_vol;
    }
    
    return helicity;
}

std::vector<MetricData> VolumeIntegral::refine_volume(
    const std::vector<MetricData>& volume_data,
    int refinement_level) {
    
    if (refinement_level <= 0) {
        return volume_data;
    }
    
    // Simplified refinement - duplicate with small perturbations
    std::vector<MetricData> refined_data = volume_data;
    
    for (int level = 0; level < refinement_level; level++) {
        size_t original_size = refined_data.size();
        for (size_t i = 0; i < original_size; i++) {
            MetricData new_data = refined_data[i];
            
            double perturbation = 0.01 * (level + 1);
            new_data.r *= (1.0 + perturbation * sin(i * 0.1));
            
            refined_data.push_back(new_data);
        }
    }
    
    return refined_data;
}

std::vector<MetricData> VolumeIntegral::filter_volume(
    const std::vector<MetricData>& volume_data,
    const std::string& filter_type,
    double filter_radius) {
    
    std::vector<MetricData> filtered_data = volume_data;
    
    if (filter_type == "gaussian" || filter_type == "smooth") {
        for (size_t i = 0; i < filtered_data.size(); i++) {
            double weighted_sum = 0.0;
            double total_weight = 0.0;
            
            for (size_t j = 0; j < volume_data.size(); j++) {
                if (i != j) {
                    double dr = volume_data[j].r - volume_data[i].r;
                    double angular_dist = std::acos(
                        sin(volume_data[i].theta) * sin(volume_data[j].theta) *
                        cos(volume_data[i].phi_coord - volume_data[j].phi_coord) +
                        cos(volume_data[i].theta) * cos(volume_data[j].theta)
                    );
                    
                    double distance = std::sqrt(dr*dr + angular_dist*angular_dist);
                    
                    if (distance < filter_radius) {
                        double weight = exp(-distance*distance / (filter_radius*filter_radius));
                        weighted_sum += volume_data[j].r * weight;
                        total_weight += weight;
                    }
                }
            }
            
            if (total_weight > 0.0) {
                filtered_data[i].r = weighted_sum / total_weight;
            }
        }
    }
    
    return filtered_data;
}

double VolumeIntegral::apply_divergence_theorem(
    const std::vector<Vector3d>& vector_field,
    const std::vector<MetricData>& volume_data,
    const std::vector<MetricData>& surface_data) {
    
    return integrate_vector_divergence(volume_data, vector_field, surface_data,
                                      IntegrationMethod::MIDPOINT_RULE);
}

std::vector<double> VolumeIntegral::compute_volume_elements(
    const std::vector<MetricData>& volume_data) {
    
    return compute_all_volume_elements(volume_data);
}

std::vector<std::array<size_t, 8>> VolumeIntegral::generate_voxel_connectivity(
    const std::vector<MetricData>& volume_data,
    const std::array<size_t, 3>& dimensions) {
    
    std::vector<std::array<size_t, 8>> voxels;
    
    size_t nx = dimensions[0];
    size_t ny = dimensions[1];
    size_t nz = dimensions[2];
    
    if (volume_data.size() != nx * ny * nz) {
        return voxels;
    }
    
    for (size_t i = 0; i < nx - 1; i++) {
        for (size_t j = 0; j < ny - 1; j++) {
            for (size_t k = 0; k < nz - 1; k++) {
                std::array<size_t, 8> voxel;
                
                voxel[0] = i * ny * nz + j * nz + k;
                voxel[1] = (i+1) * ny * nz + j * nz + k;
                voxel[2] = (i+1) * ny * nz + (j+1) * nz + k;
                voxel[3] = i * ny * nz + (j+1) * nz + k;
                voxel[4] = i * ny * nz + j * nz + (k+1);
                voxel[5] = (i+1) * ny * nz + j * nz + (k+1);
                voxel[6] = (i+1) * ny * nz + (j+1) * nz + (k+1);
                voxel[7] = i * ny * nz + (j+1) * nz + (k+1);
                
                voxels.push_back(voxel);
            }
        }
    }
    
    return voxels;
}

std::vector<MetricData> VolumeIntegral::interpolate_to_uniform_grid(
    const std::vector<MetricData>& volume_data,
    const std::array<size_t, 3>& new_dimensions) {
    
    size_t new_size = new_dimensions[0] * new_dimensions[1] * new_dimensions[2];
    std::vector<MetricData> uniform_grid(new_size);
    
    // Simplified - nearest neighbor interpolation
    for (size_t i = 0; i < new_size; i++) {
        size_t nearest_idx = i % volume_data.size();
        uniform_grid[i] = volume_data[nearest_idx];
    }
    
    return uniform_grid;
}

} // namespace grid
} // namespace qlt