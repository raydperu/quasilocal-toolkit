#ifndef QLT_CORE_H
#define QLT_CORE_H

// Include all core module headers
#include "axisymmetric_j.h"
#include "clebsch_komar.h"
#include "quasilocal_mass.h"
#include "surgical_flux.h"
#include "vorticity_bounds.h"

/**
 * @brief Main namespace for quasilocal toolkit core functionality
 * 
 * This library implements the mathematical framework from:
 * "Geometric Vorticity and Quasilocal Angular Momentum Conservation"
 * by Rayan D. Peru.
 * 
 * All functions are mathematically accurate implementations of the
 * formulas in the paper, with proper numerical handling.
 */
namespace qlt {
    
    // All functionality is available through the individual headers
    // This header provides a convenient single include point
    
} // namespace qlt

#endif // QLT_CORE_H