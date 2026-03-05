#ifndef QLT_VOLUME_INTEGRAL_H
#define QLT_VOLUME_INTEGRAL_H

#include "../core/quasilocal_quantities.h"
#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <functional>
#include <string>

namespace qlt {
namespace grid {

/**
 * @brief Volume integration utilities
 */
class VolumeIntegral {
public:
    /**
     * @brief Integration methods for 3D volumes
     */
    enum class IntegrationMethod {
        MIDPOINT_RULE,      // Simple midpoint rule
        TRAPEZOIDAL_RULE,   // Trapezoidal rule
        SIMPSONS_RULE,      // Simpson's rule
        MONTE_CARLO,        // Monte Carlo integration
        ADAPTIVE_CUBATURE   // Adaptive cubature
    };
    
    /**
     * @brief Compute volume integral of scalar function
     */
    static double integrate_scalar(
        const std::vector<MetricData>& volume_data,
        const std::vector<double>& scalar_values,
        IntegrationMethod method = IntegrationMethod::MIDPOINT_RULE);
    
    /**
     * @brief Compute volume integral of vector field (divergence theorem)
     */
    static double integrate_vector_divergence(
        const std::vector<MetricData>& volume_data,
        const std::vector<Vector3d>& vector_field,
        const std::vector<MetricData>& surface_data,
        IntegrationMethod method = IntegrationMethod::MIDPOINT_RULE);
    
    /**
     * @brief Compute Clebsch-Komar volume term
     *        ∫_Σ [α dβ ∧ dξ^♭ + Φ d(dβ ∧ ξ^♭)]
     */
    static double compute_clebsch_komar_volume_integral(
        const std::vector<MetricData>& volume_data);
    
    /**
     * @brief Compute volume of region
     */
    static double compute_volume(
        const std::vector<MetricData>& volume_data);
    
    /**
     * @brief Compute center of mass
     */
    static std::array<double, 3> compute_center_of_mass(
        const std::vector<MetricData>& volume_data,
        const std::vector<double>& density = {});
    
    /**
     * @brief Extract subvolume from larger grid
     */
    static std::vector<MetricData> extract_subvolume(
        const std::vector<std::vector<std::vector<MetricData>>>& grid_data,
        const std::array<double, 3>& min_corner,
        const std::array<double, 3>& max_corner);
    
    /**
     * @brief Compute integral for vorticity bound
     *        ∫_Σ ||α||_g ||dβ||_g² √h d³x
     */
    static double compute_vorticity_bound_integral(
        const std::vector<MetricData>& volume_data);
    
    /**
     * @brief Compute helicity integral ∫ V ∧ dV
     */
    static double compute_helicity_integral(
        const std::vector<MetricData>& volume_data);
    
    /**
     * @brief Refine volume grid (increase resolution)
     */
    static std::vector<MetricData> refine_volume(
        const std::vector<MetricData>& volume_data,
        int refinement_level = 1);
    
    /**
     * @brief Apply filter to volume data (smoothing)
     */
    static std::vector<MetricData> filter_volume(
        const std::vector<MetricData>& volume_data,
        const std::string& filter_type = "gaussian",
        double filter_radius = 1.0);
    
    /**
     * @brief Compute volume integral using divergence theorem
     *        ∫_V ∇·F dV = ∫_∂V F·n dA
     */
    static double apply_divergence_theorem(
        const std::vector<Vector3d>& vector_field,
        const std::vector<MetricData>& volume_data,
        const std::vector<MetricData>& surface_data);
    
private:
    /**
     * @brief Compute volume element dV for each point
     */
    static std::vector<double> compute_volume_elements(
        const std::vector<MetricData>& volume_data);
    
    /**
     * @brief Generate 3D grid connectivity
     */
    static std::vector<std::array<size_t, 8>> generate_voxel_connectivity(
        const std::vector<MetricData>& volume_data,
        const std::array<size_t, 3>& dimensions);
    
    /**
     * @brief Interpolate field to uniform grid
     */
    static std::vector<MetricData> interpolate_to_uniform_grid(
        const std::vector<MetricData>& volume_data,
        const std::array<size_t, 3>& new_dimensions);
};

} // namespace grid
} // namespace qlt

#endif // QLT_VOLUME_INTEGRAL_H