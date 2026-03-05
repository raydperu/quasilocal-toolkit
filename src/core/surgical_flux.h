#ifndef QLT_SURGICAL_FLUX_H
#define QLT_SURGICAL_FLUX_H

#include "quasilocal_quantities.h"

namespace qlt {

/**
 * @brief Compute surgical flux (angular momentum change across excision)
 */
double compute_surgical_flux(
    const std::vector<MetricData>& inner_data,
    const std::vector<MetricData>& outer_data,
    bool axisymmetric = false);

/**
 * @brief Compute cumulative angular momentum loss through multiple surgeries
 */
double compute_cumulative_angular_momentum_loss(
    const std::vector<std::vector<MetricData>>& surgical_boundaries,
    const std::vector<bool>& axisymmetric_flags);

/**
 * @brief Verify angular momentum conservation across surgery
 */
bool verify_surgery_conservation(
    double j_before,
    double j_after,
    const std::vector<MetricData>& inner_boundary,
    const std::vector<MetricData>& outer_boundary,
    double tolerance = 1e-6);

} // namespace qlt

#endif // QLT_SURGICAL_FLUX_H
