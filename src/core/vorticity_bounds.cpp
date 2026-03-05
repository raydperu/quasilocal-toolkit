#include "quasilocal_quantities.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace qlt {

namespace { // Helper functions

/**
 * @brief Compute vorticity 1-form V = dΦ + α dβ
 */
Vector3d compute_vorticity_1form(const MetricData& data) {
    return data.dphi_dx + data.alpha * data.dbeta_dx;
}

/**
 * @brief Compute exterior derivative of vorticity: dV
 *        Using Clebsch constraint: dV = dα ∧ dβ
 */
Matrix3d compute_dvorticity(const MetricData& data) {
    // dV as antisymmetric 2-form
    // dV_{ij} = ∂_i V_j - ∂_j V_i = (∂_i α)(∂_j β) - (∂_j α)(∂_i β)
    
    Matrix3d dv;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            dv(i, j) = data.dalpha_dx[i] * data.dbeta_dx[j] - 
                       data.dalpha_dx[j] * data.dbeta_dx[i];
        }
    }
    
    return dv;
}

/**
 * @brief Compute norm of 1-form with respect to metric: ||ω||_g = √(h^{ij} ω_i ω_j)
 */
double compute_form_norm(const Vector3d& form, const Matrix3d& h_ij) {
    Matrix3d h_inv = h_ij.inverse();
    double norm_sq = form.transpose() * h_inv * form;
    return std::sqrt(std::max(0.0, norm_sq));
}

/**
 * @brief Compute norm of 2-form: ||Ω||_g
 */
double compute_2form_norm(const Matrix3d& form_2, const Matrix3d& h_ij) {
    // For antisymmetric 2-form Ω_{ij}
    // ||Ω||² = (1/2) h^{ik} h^{jl} Ω_{ij} Ω_{kl}
    
    Matrix3d h_inv = h_ij.inverse();
    double norm_sq = 0.0;
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                for (int l = 0; l < 3; ++l) {
                    norm_sq += h_inv(i, k) * h_inv(j, l) * 
                               form_2(i, j) * form_2(k, l);
                }
            }
        }
    }
    
    return std::sqrt(std::max(0.0, norm_sq / 2.0));
}

} // anonymous namespace

std::pair<double, double> compute_vorticity_bound(const MetricData& data) {
    // Theorem 3.1/4.1: ||V ∧ dV||_g ≤ ||α||_g ||dβ||_g²
    
    // Compute left-hand side: ||V ∧ dV||_g
    Vector3d V = compute_vorticity_1form(data);
    Matrix3d dV = compute_dvorticity(data);
    
    // V ∧ dV is a 3-form in 3D (fully antisymmetric)
    // In 3D, wedge product of 1-form and 2-form gives a 3-form
    // The norm of a 3-form in 3D is the absolute value of its component
    
    // Compute component: ε^{ijk} V_i dV_{jk} / 3!
    // In 3D, this is essentially the triple product
    
    double v_wedge_dv_component = 
        V[0] * (dV(1, 2) - dV(2, 1)) +
        V[1] * (dV(2, 0) - dV(0, 2)) +
        V[2] * (dV(0, 1) - dV(1, 0));
    
    v_wedge_dv_component /= 6.0; // 3! = 6
    
    // Norm of 3-form in 3D: ||ω|| = |ω_{123}| / √|det h|
    double lhs = std::abs(v_wedge_dv_component) / std::sqrt(data.h_ij.determinant());
    
    // Compute right-hand side: ||α||_g ||dβ||_g²
    // ||α||_g = |α| (scalar norm)
    double alpha_norm = std::abs(data.alpha);
    
    // ||dβ||_g = norm of 1-form dβ
    double dbeta_norm = compute_form_norm(data.dbeta_dx, data.h_ij);
    
    double rhs = alpha_norm * dbeta_norm * dbeta_norm;
    
    return std::make_pair(lhs, rhs);
}

/**
 * @brief Check if vorticity bound is satisfied
 */
bool check_vorticity_bound(const MetricData& data, double tolerance = 1e-10) {
    auto [lhs, rhs] = compute_vorticity_bound(data);
    bool satisfied = (lhs <= rhs + tolerance);
    
    if (!satisfied) {
        std::cerr << "[WARNING] Vorticity bound violated: " 
                  << lhs << " > " << rhs << std::endl;
    }
    
    return satisfied;
}

/**
 * @brief Compute integrated angular momentum bound (Theorem 3.2/4.2)
 *        |J[ξ]| ≤ (1/16π) ∫_Σ ||α||_g ||dβ||_g² √h d³x + C(M_ADM)
 */
std::pair<double, double> compute_integrated_angular_momentum_bound(
    const std::vector<MetricData>& volume_data,
    double adm_mass,
    double mass_constant_factor = 1.0) {
    
    double volume_integral = 0.0;
    
    for (const auto& data : volume_data) {
        // Compute integrand: ||α||_g ||dβ||_g² √h
        double alpha_norm = std::abs(data.alpha);
        double dbeta_norm = compute_form_norm(data.dbeta_dx, data.h_ij);
        double sqrt_h = std::sqrt(data.h_ij.determinant());
        
        double integrand = alpha_norm * dbeta_norm * dbeta_norm * sqrt_h;
        volume_integral += integrand;
    }
    
    // First term: (1/16π) times volume integral
    double first_term = volume_integral / (16.0 * M_PI);
    
    // Second term: C(M_ADM) - simplified as constant times mass
    // In paper: C(M_ADM) depends on ADM mass
    double second_term = mass_constant_factor * adm_mass;
    
    double bound = first_term + second_term;
    
    return std::make_pair(first_term, bound);
}

/**
 * @brief Compute helicity ∫ V ∧ dV (Theorem 2.4)
 */
double compute_helicity(const std::vector<MetricData>& volume_data) {
    double helicity = 0.0;
    
    for (const auto& data : volume_data) {
        Vector3d V = compute_vorticity_1form(data);
        Matrix3d dV = compute_dvorticity(data);
        
        // V ∧ dV component (3-form)
        double component = 
            V[0] * (dV(1, 2) - dV(2, 1)) +
            V[1] * (dV(2, 0) - dV(0, 2)) +
            V[2] * (dV(0, 1) - dV(1, 0));
        
        component /= 6.0; // 3! = 6
        
        // Volume element √h d³x
        double dV_volume = std::sqrt(data.h_ij.determinant());
        
        helicity += component * dV_volume;
    }
    
    return helicity;
}

/**
 * @brief Check zero angular momentum condition (Corollary 3.3/4.3)
 */
std::pair<bool, bool> check_zero_angular_momentum(
    const std::vector<MetricData>& data,
    double tolerance) {
    
    bool alpha_zero = true;
    bool dbeta_zero = true;
    
    double max_alpha = 0.0;
    double max_dbeta_norm = 0.0;
    
    for (const auto& d : data) {
        max_alpha = std::max(max_alpha, std::abs(d.alpha));
        double dbeta_norm = compute_form_norm(d.dbeta_dx, d.h_ij);
        max_dbeta_norm = std::max(max_dbeta_norm, dbeta_norm);
    }
    
    alpha_zero = (max_alpha < tolerance);
    dbeta_zero = (max_dbeta_norm < tolerance);
    
    std::cout << "[DEBUG] Zero angular momentum check:" << std::endl;
    std::cout << "  max|α| = " << max_alpha 
              << " (zero if < " << tolerance << ")" << std::endl;
    std::cout << "  max||dβ|| = " << max_dbeta_norm
              << " (zero if < " << tolerance << ")" << std::endl;
    std::cout << "  α ≡ 0? " << (alpha_zero ? "YES" : "NO") << std::endl;
    std::cout << "  dβ ≡ 0? " << (dbeta_zero ? "YES" : "NO") << std::endl;
    
    return std::make_pair(alpha_zero, dbeta_zero);
}

} // namespace qlt
