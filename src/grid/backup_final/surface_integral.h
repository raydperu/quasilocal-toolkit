#ifndef QLT_SURFACE_INTEGRAL_H
#define QLT_SURFACE_INTEGRAL_H

#include "../core/quasilocal_quantities.h"
#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <functional>

namespace qlt {
namespace grid {

/**
 * @brief Surface geometry and integration utilities
 */
class SurfaceIntegral {
public:
    /**
     * @brief Integration methods for 2D surfaces
     */
    enum class IntegrationMethod {
        MIDPOINT_RULE,      // Simple midpoint rule
        TRAPEZOIDAL_RULE,   // Trapezoidal rule
        SIMPSONS_RULE,      // Simpson's rule (for smooth functions)
        SPHERICAL_QUADRATURE, // Specialized for spheres
        ADAPTIVE_QUADRATURE   // Adaptive refinement
    };
    
    /**
     * @brief Compute surface integral of scalar function
     */
    static double integrate_scalar(
        const std::vector<MetricData>& surface_data,
        const std::vector<double>& scalar_values,
        IntegrationMethod method = IntegrationMethod::MIDPOINT_RULE);
    
    /**
     * @brief Compute surface integral of vector field (flux)
     */
    static double integrate_vector_flux(
        const std::vector<MetricData>& surface_data,
        const std::vector<Vector3d>& vector_field,
        IntegrationMethod method = IntegrationMethod::MIDPOINT_RULE);
    
    /**
     * @brief Compute Komar surface integral: ∫_S ∇^μ ξ^ν dS_{μν}
     */
    static double compute_komar_surface_integral(
        const std::vector<MetricData>& surface_data);
    
    /**
     * @brief Compute area of 2-surface
     */
    static double compute_surface_area(
        const std::vector<MetricData>& surface_data);
    
    /**
     * @brief Compute centroid of surface
     */
    static std::array<double, 3> compute_surface_centroid(
        const std::vector<MetricData>& surface_data);
    
    /**
     * @brief Compute surface normal vectors (if not provided)
     */
    static std::vector<Vector3d> compute_surface_normals(
        const std::vector<MetricData>& surface_data);
    
    /**
     * @brief Compute surface curvature invariants
     */
    static std::vector<double> compute_surface_curvature(
        const std::vector<MetricData>& surface_data);
    
    /**
     * @brief Check if surface is closed (no boundary)
     */
    static bool is_closed_surface(
        const std::vector<MetricData>& surface_data);
    
    /**
     * @brief Extract surface from 3D grid at specified radius
     */
    static std::vector<MetricData> extract_spherical_surface(
        const std::vector<std::vector<std::vector<MetricData>>>& grid_data,
        double radius,
        size_t num_theta = 64,
        size_t num_phi = 128);
    
    /**
     * @brief Extract surface from 3D grid using isosurface extraction
     */
    static std::vector<MetricData> extract_isosurface(
        const std::vector<std::vector<std::vector<MetricData>>>& grid_data,
        const std::string& field_name,  // e.g., "r", "density"
        double isovalue,
        bool use_marching_cubes = true);
    
    /**
     * @brief Refine surface mesh (increase resolution)
     */
    static std::vector<MetricData> refine_surface(
        const std::vector<MetricData>& surface_data,
        int refinement_level = 1);
    
    /**
     * @brief Smooth surface (reduce noise)
     */
    static std::vector<MetricData> smooth_surface(
        const std::vector<MetricData>& surface_data,
        int smoothing_iterations = 3,
        double smoothing_factor = 0.5);
    
    /**
     * @brief Compute surface integral for surgical flux
     *        ∫_∂V α dβ ∧ ξ^♭
     */
    static double compute_surgical_flux_integral(
        const std::vector<MetricData>& surface_data);
    
private:
    /**
     * @brief Compute surface area element dA for each point
     */
    static std::vector<double> compute_area_elements(
        const std::vector<MetricData>& surface_data);
    
    /**
     * @brief Generate connectivity information for surface mesh
     */
    static std::vector<std::array<size_t, 3>> generate_face_connectivity(
        const std::vector<MetricData>& surface_data);
    
    /**
     * @brief Interpolate missing normal vectors
     */
    static void interpolate_normals(
        std::vector<MetricData>& surface_data);
};

} // namespace grid
} // namespace qlt

#endif // QLT_SURFACE_INTEGRAL_H