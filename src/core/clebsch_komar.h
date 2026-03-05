#ifndef QLT_CLEBSCH_KOMAR_H
#define QLT_CLEBSCH_KOMAR_H

#include "qlt/quasilocal_quantities.h"

namespace qlt {

/**
 * @brief Core Clebsch-Komar angular momentum computation
 */
double compute_clebsch_komar_j(
    const std::vector<MetricData>& surface_data,
    const std::vector<MetricData>& volume_data,
    double kappa = 2.0);

/**
 * @brief Traditional Komar surface integral (for comparison)
 */
double compute_komar_surface_integral(
    const std::vector<MetricData>& surface_data);

/**
 * @brief Verify Clebsch constraint C_{μν} = 0
 */
bool verify_clebsch_constraint(const MetricData& data, double tolerance = 1e-8);

} // namespace qlt

#endif // QLT_CLEBSCH_KOMAR_H