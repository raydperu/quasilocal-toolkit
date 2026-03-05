#ifndef QLT_KILLING_VECTOR_H
#define QLT_KILLING_VECTOR_H

#include "qlt/quasilocal_quantities.h"
#include <vector>
#include <string>
#include <memory>

namespace qlt {
namespace grid {

/**
 * @brief Killing vector field types
 */
enum class KillingVectorType {
    ROTATIONAL_Z,      // Rotation about z-axis: ξ = ∂/∂φ
    ROTATIONAL_X,      // Rotation about x-axis
    ROTATIONAL_Y,      // Rotation about y-axis
    TRANSLATIONAL_X,   // Translation along x: ξ = ∂/∂x
    TRANSLATIONAL_Y,   // Translation along y
    TRANSLATIONAL_Z,   // Translation along z
    BOOST_X,          // Lorentz boost along x
    BOOST_Y,          // Lorentz boost along y
    BOOST_Z,          // Lorentz boost along z
    CUSTOM,           // User-defined Killing vector
    APPROXIMATE       // Approximate Killing vector from metric
};

/**
 * @brief Class for generating and working with Killing vector fields
 */
class KillingVectorGenerator {
public:
    KillingVectorGenerator();
    
    /**
     * @brief Generate Killing vector field of specified type
     */
    std::vector<MetricData> generate_killing_field(
        const std::vector<MetricData>& grid_data,
        KillingVectorType type,
        const std::array<double, 3>& parameters = {0.0, 0.0, 0.0});
    
    /**
     * @brief Generate rotational Killing vector about specified axis
     */
    std::vector<MetricData> generate_rotational_killing_vector(
        const std::vector<MetricData>& grid_data,
        const std::array<double, 3>& axis = {0.0, 0.0, 1.0},  // Default: z-axis
        double angular_velocity = 1.0);
    
    /**
     * @brief Generate translational Killing vector along specified direction
     */
    std::vector<MetricData> generate_translational_killing_vector(
        const std::vector<MetricData>& grid_data,
        const std::array<double, 3>& direction = {1.0, 0.0, 0.0});
    
    /**
     * @brief Compute approximate Killing vector from metric using eigen-decomposition
     *        Finds vector field that approximately satisfies Killing equation
     */
    std::vector<MetricData> compute_approximate_killing_vector(
        const std::vector<MetricData>& grid_data,
        int max_iterations = 100,
        double tolerance = 1e-6);
    
    /**
     * @brief Check if vector field satisfies Killing equation
     *        ∇_i ξ_j + ∇_j ξ_i = 0
     */
    std::pair<double, double> check_killing_equation(
        const std::vector<MetricData>& data_with_xi,
        double& max_violation,
        double& rms_violation);
    
    /**
     * @brief Compute Lie derivative of metric along vector field
     */
    std::vector<Matrix3d> compute_lie_derivative(
        const std::vector<MetricData>& data);
    
    /**
     * @brief Generate Killing vectors for common spacetimes
     */
    std::vector<MetricData> generate_for_spacetime(
        const std::vector<MetricData>& grid_data,
        const std::string& spacetime_name,  // "Kerr", "Schwarzschild", "Minkowski"
        const std::map<std::string, double>& parameters);
    
    /**
     * @brief Read Killing vector from file
     */
    std::vector<MetricData> read_from_file(
        const std::string& filename,
        const std::vector<MetricData>& grid_data);
    
    /**
     * @brief Save Killing vector to file
     */
    void save_to_file(
        const std::vector<MetricData>& data_with_xi,
        const std::string& filename);
    
private:
    /**
     * @brief Compute derivatives of Killing vector field using finite differences
     */
    Matrix3d compute_killing_derivatives(
        const MetricData& data,
        const std::vector<MetricData>& neighborhood);
    
    /**
     * @brief Find Killing vector that minimizes Killing equation violation
     */
    Vector3d optimize_killing_vector(
        const MetricData& data,
        const std::vector<MetricData>& neighborhood);
};

} // namespace grid
} // namespace qlt

#endif // QLT_KILLING_VECTOR_H