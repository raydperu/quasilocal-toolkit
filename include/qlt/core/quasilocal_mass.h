#ifndef QLT_QUASILOCAL_MASS_H
#define QLT_QUASILOCAL_MASS_H

#include "qlt/quasilocal_quantities.h"

namespace qlt {

/**
 * @brief Compute quasilocal mass via modified Brown-York formalism
 */
double compute_quasilocal_mass(
    const std::vector<MetricData>& surface_data,
    bool use_reference = true);

/**
 * @brief Approximate ADM mass from asymptotic data
 */
double compute_adm_mass_approximation(
    const std::vector<MetricData>& surface_data);

/**
 * @brief Compute mass aspect function M(r) for spherical symmetry
 */
std::vector<double> compute_mass_aspect(
    const std::vector<MetricData>& surface_data_at_radii);

} // namespace qlt

#endif // QLT_QUASILOCAL_MASS_H