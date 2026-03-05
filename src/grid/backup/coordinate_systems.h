#ifndef QLT_COORDINATE_SYSTEMS_H
#define QLT_COORDINATE_SYSTEMS_H

#include "qlt/quasilocal_quantities.h"
#include <vector>
#include <array>
#include <cmath>

namespace qlt {
namespace grid {

/**
 * @brief Coordinate transformation utilities
 */
class CoordinateSystems {
public:
    /**
     * @brief Convert Cartesian to spherical coordinates
     */
    static std::array<double, 3> cartesian_to_spherical(double x, double y, double z);
    
    /**
     * @brief Convert spherical to Cartesian coordinates
     */
    static std::array<double, 3> spherical_to_cartesian(double r, double theta, double phi);
    
    /**
     * @brief Convert Cartesian to cylindrical coordinates
     */
    static std::array<double, 3> cartesian_to_cylindrical(double x, double y, double z);
    
    /**
     * @brief Convert cylindrical to Cartesian coordinates
     */
    static std::array<double, 3> cylindrical_to_cartesian(double rho, double phi, double z);
    
    /**
     * @brief Compute Jacobian matrix for coordinate transformation
     * @param from Coordinate system to transform from
     * @param to Coordinate system to transform to
     * @param point Coordinates in 'from' system
     * @return 3x3 Jacobian matrix ∂(to)/∂(from)
     */
    static Matrix3d compute_jacobian(
        const std::string& from,
        const std::string& to,
        const std::array<double, 3>& point);
    
    /**
     * @brief Transform vector from one coordinate system to another
     */
    static Vector3d transform_vector(
        const Vector3d& vector,
        const std::string& from,
        const std::string& to,
        const std::array<double, 3>& point);
    
    /**
     * @brief Transform metric from one coordinate system to another
     */
    static Matrix3d transform_metric(
        const Matrix3d& metric,
        const std::string& from,
        const std::string& to,
        const std::array<double, 3>& point);
    
    /**
     * @brief Generate spherical grid coordinates
     * @param r_min Minimum radius
     * @param r_max Maximum radius
     * @param num_r Number of radial points
     * @param num_theta Number of θ points
     * @param num_phi Number of φ points
     * @return Grid of (r, θ, φ) coordinates
     */
    static std::vector<std::vector<std::vector<std::array<double, 3>>>>
    generate_spherical_grid(
        double r_min, double r_max, size_t num_r,
        size_t num_theta, size_t num_phi);
    
    /**
     * @brief Generate Cartesian grid coordinates
     */
    static std::vector<std::vector<std::vector<std::array<double, 3>>>>
    generate_cartesian_grid(
        double x_min, double x_max, size_t num_x,
        double y_min, double y_max, size_t num_y,
        double z_min, double z_max, size_t num_z);
    
    /**
     * @brief Interpolate data from one grid to another
     */
    static std::vector<MetricData> interpolate_to_grid(
        const std::vector<MetricData>& source_data,
        const std::vector<std::array<double, 3>>& target_coords,
        const std::string& interpolation_method = "trilinear");
    
    /**
     * @brief Check if coordinates are valid (e.g., r >= 0, 0 <= θ <= π)
     */
    static bool validate_coordinates(
        const std::array<double, 3>& coords,
        const std::string& system);
    
    /**
     * @brief Compute coordinate basis vectors
     * @param system Coordinate system
     * @param point Coordinate point
     * @return Array of 3 basis vectors
     */
    static std::array<Vector3d, 3> compute_basis_vectors(
        const std::string& system,
        const std::array<double, 3>& point);
    
    /**
     * @brief Compute Christoffel symbols for a metric in given coordinates
     */
    static std::array<std::array<std::array<double, 3>, 3>, 3>
    compute_christoffel_symbols(
        const Matrix3d& metric,
        const std::string& coordinate_system,
        const std::array<double, 3>& point,
        double epsilon = 1e-6);
    
    /**
     * @brief Compute covariant derivative of a vector
     */
    static Matrix3d compute_covariant_derivative(
        const Vector3d& vector,
        const Matrix3d& metric,
        const std::string& coordinate_system,
        const std::array<double, 3>& point,
        double epsilon = 1e-6);
};

} // namespace grid
} // namespace qlt

#endif // QLT_COORDINATE_SYSTEMS_H