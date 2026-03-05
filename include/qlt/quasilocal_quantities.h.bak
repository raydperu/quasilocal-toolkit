#ifndef QLT_QUASILOCAL_QUANTITIES_H
#define QLT_QUASILOCAL_QUANTITIES_H

#include <vector>
#include <array>
#include <Eigen/Dense>

namespace qlt {

/**
 * @brief Core quasilocal quantity computations based on
 * "Geometric Vorticity and Quasilocal Angular Momentum Conservation"
 * by Rayan D. Peru.
 */

// ============================================================================
// TYPE ALIASES
// ============================================================================

using Vector3d = Eigen::Vector3d;
using Matrix3d = Eigen::Matrix3d;

/**
 * @brief Structure holding metric data at a point
 */
struct MetricData {
    Matrix3d h_ij;          // Spatial metric (3x3 symmetric positive definite)
    Matrix3d K_ij;          // Extrinsic curvature (3x3 symmetric)
    Vector3d xi;            // Killing vector field components
    Matrix3d dxi_dx;        // Derivative ∂_i ξ_j (3x3)
    double sqrt_det_h;      // √det(h_ij)
    
    // Clebsch potentials and derivatives
    double alpha;           // Clebsch potential α
    double phi;             // Clebsch potential Φ
    double beta;            // Clebsch potential β
    Vector3d dalpha_dx;     // ∇α
    Vector3d dbeta_dx;      // ∇β
    Vector3d dphi_dx;       // ∇Φ
    
    // Normal vector to 2-surface (if on boundary)
    Vector3d normal;        // Unit normal (outward pointing)
    
    // Position information
    double r;               // Radial coordinate (if axisymmetric)
    double theta;           // Polar angle
    double phi_coord;       // Azimuthal angle
};

// ============================================================================
// CORE FUNCTIONS
// ============================================================================

/**
 * @brief Compute quasilocal angular momentum in axisymmetric case
 *        Theorem 4.3: J(r) = (1/3) α(r) r² f(r)
 *        where β = φ + f(r)cosθ
 * 
 * @param alpha α(r) values
 * @param f f(r) = β - φ (phase function)
 * @param r Radial coordinates
 * @return std::vector<double> J(r) values
 */
std::vector<double> compute_axisymmetric_j(
    const std::vector<double>& alpha,
    const std::vector<double>& f,
    const std::vector<double>& r);

/**
 * @brief Compute Clebsch-Komar angular momentum for arbitrary 2-surface S
 *        Eq. 1.3: J_S[ξ] = (1/16π) ∫_S ∇^μ ξ^ν dS_μν 
 *                       + (κ/16π) ∫_Σ_S [α dβ ∧ dξ^♭ + Φ d(dβ ∧ ξ^♭)]
 * 
 * @param surface_data Data on 2-surface S (boundary)
 * @param volume_data Data in volume Σ_S (inside S)
 * @param kappa Coupling constant κ (default = 2 as in paper)
 * @return double Angular momentum J_S[ξ]
 */
double compute_clebsch_komar_j(
    const std::vector<MetricData>& surface_data,
    const std::vector<MetricData>& volume_data,
    double kappa = 2.0);

/**
 * @brief Compute quasilocal mass via modified Brown-York formalism
 *        Adapted to be consistent with Clebsch framework
 * 
 * @param surface_data Data on 2-surface S
 * @param use_reference Use flat space reference curvature k₀
 * @return double Quasilocal mass M_S
 */
double compute_quasilocal_mass(
    const std::vector<MetricData>& surface_data,
    bool use_reference = true);

/**
 * @brief Compute surgical flux across excised region
 *        Theorem 4.4: ΔJ = J(λ⁺) - J(λ⁻) = -(1/8π) ∫_∂V α dβ ∧ ξ^♭
 *        Corollary 4.6 (axisymmetric): ΔJ = (1/3)[α(r₂)r₂²f(r₂) - α(r₁)r₁²f(r₁)]
 * 
 * @param inner_data Data on inner boundary S₁
 * @param outer_data Data on outer boundary S₂
 * @param axisymmetric Use axisymmetric simplification
 * @return double Angular momentum change ΔJ
 */
double compute_surgical_flux(
    const std::vector<MetricData>& inner_data,
    const std::vector<MetricData>& outer_data,
    bool axisymmetric = false);

/**
 * @brief Check if zero angular momentum condition is satisfied
 *        Corollary 3.3/4.3: If α ≡ 0 or dβ ≡ 0, then J = 0
 * 
 * @param data Metric data
 * @param tolerance Numerical tolerance
 * @return std::pair<bool, bool> (alpha_zero, dbeta_zero)
 */
std::pair<bool, bool> check_zero_angular_momentum(
    const std::vector<MetricData>& data,
    double tolerance = 1e-10);

/**
 * @brief Compute vorticity bound from Theorem 3.1/4.1
 *        ||V ∧ dV||_g ≤ ||α||_g ||dβ||_g²
 *        Returns (lhs, rhs) for verification
 */
std::pair<double, double> compute_vorticity_bound(
    const MetricData& data);

} // namespace qlt

#endif // QLT_QUASILOCAL_QUANTITIES_H