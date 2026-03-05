#ifndef QLT_VORTICITY_BOUNDS_H
#define QLT_VORTICITY_BOUNDS_H

#include "quasilocal_quantities.h"

namespace qlt {

/**
 * @brief Compute vorticity bound from Theorem 3.1/4.1
 *        Returns (||V ∧ dV||_g, ||α||_g ||dβ||_g²)
 */
std::pair<double, double> compute_vorticity_bound(const MetricData& data);

/**
 * @brief Check if vorticity bound is satisfied
 */
bool check_vorticity_bound(const MetricData& data, double tolerance = 1e-10);

/**
 * @brief Compute integrated angular momentum bound (Theorem 3.2/4.2)
 */
std::pair<double, double> compute_integrated_angular_momentum_bound(
    const std::vector<MetricData>& volume_data,
    double adm_mass,
    double mass_constant_factor = 1.0);

/**
 * @brief Compute helicity ∫ V ∧ dV
 */
double compute_helicity(const std::vector<MetricData>& volume_data);

/**
 * @brief Check zero angular momentum condition (Corollary 3.3/4.3)
 */
std::pair<bool, bool> check_zero_angular_momentum(
    const std::vector<MetricData>& data,
    double tolerance = 1e-10);

} // namespace qlt

#endif // QLT_VORTICITY_BOUNDS_H
