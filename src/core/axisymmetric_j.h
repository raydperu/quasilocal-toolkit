#ifndef QLT_AXISYMMETRIC_J_H
#define QLT_AXISYMMETRIC_J_H

#include "quasilocal_quantities.h"

namespace qlt {

// Core functions from axisymmetric_j.cpp
std::vector<double> compute_axisymmetric_j(
    const std::vector<double>& alpha,
    const std::vector<double>& f,
    const std::vector<double>& r);

double compute_axisymmetric_surgical_flux(
    double alpha1, double r1, double f1,
    double alpha2, double r2, double f2);

// Helper functions for extracting from MetricData
std::vector<double> extract_alpha_from_surface(
    const std::vector<MetricData>& surface_data);

std::vector<double> extract_f_from_surface(
    const std::vector<MetricData>& surface_data);

std::vector<double> extract_r_from_surface(
    const std::vector<MetricData>& surface_data);

} // namespace qlt

#endif // QLT_AXISYMMETRIC_J_H
