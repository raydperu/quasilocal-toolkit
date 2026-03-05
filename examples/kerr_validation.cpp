/**
 * @example kerr_validation.cpp
 * 
 * Validate quasilocal angular momentum computation against Kerr metric.
 * 
 * This example:
 * 1. Creates Kerr metric data at various radii
 * 2. Computes angular momentum using axisymmetric formula
 * 3. Compares with expected value J = M*a
 * 4. Demonstrates vorticity bounds
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "qlt/quasilocal_quantities.h"
#include "grid/coordinate_systems.h"

using namespace qlt;

/**
 * @brief Create Kerr metric data at given radius
 */
MetricData create_kerr_data_at_radius(double r, double theta, 
                                     double M, double a) {
    MetricData data;
    
    data.r = r;
    data.theta = theta;
    data.phi_coord = 0.0;  // φ = 0 for simplicity
    
    // Kerr metric in Boyer-Lindquist coordinates
    // (simplified for example - full metric would be more complex)
    double Sigma = r*r + a*a*cos(theta)*cos(theta);
    double Delta = r*r - 2.0*M*r + a*a;
    
    // Set Clebsch potentials for Kerr
    data.alpha = -2.0 * M * a * r * sin(theta)*sin(theta) / Sigma;
    data.phi = 0.0;  // Constant potential
    data.beta = data.phi_coord;  // β = φ for axisymmetry
    
    // Simple flat metric for demonstration
    // (In real application, would use actual Kerr metric)
    data.h_ij = Matrix3d::Identity();
    data.sqrt_det_h = 1.0;
    
    // Rotational Killing vector about z-axis
    double x = r * sin(theta) * cos(data.phi_coord);
    double y = r * sin(theta) * sin(data.phi_coord);
    data.xi = Vector3d(-y, x, 0.0);
    
    // Derivatives (simplified)
    data.dxi_dx = Matrix3d::Zero();
    data.dxi_dx(0,1) = -1.0;
    data.dxi_dx(1,0) = 1.0;
    
    // Clebsch derivatives
    data.dalpha_dx = Vector3d::Zero();
    data.dbeta_dx = Vector3d(0.0, 0.0, 1.0);  // ∇β = ∇φ
    data.dphi_dx = Vector3d::Zero();
    
    return data;
}

int main() {
    std::cout << "=============================================\n";
    std::cout << "Kerr Metric Validation Example\n";
    std::cout << "=============================================\n\n";
    
    // Parameters
    const double M = 1.0;      // Mass
    const double a = 0.5;      // Spin parameter
    const double J_expected = M * a;
    
    std::cout << "Parameters:\n";
    std::cout << "  Mass M = " << M << "\n";
    std::cout << "  Spin a = " << a << "\n";
    std::cout << "  Expected angular momentum J = M*a = " 
              << J_expected << "\n\n";
    
    // Test at various radii
    std::vector<double> radii = {10.0, 20.0, 50.0, 100.0, 200.0};
    const double theta = M_PI / 2.0;  // Equatorial plane
    
    std::cout << std::setw(10) << "Radius" 
              << std::setw(15) << "α(r)" 
              << std::setw(15) << "f(r)" 
              << std::setw(15) << "J(r)" 
              << std::setw(15) << "Error" 
              << std::setw(15) << "Rel. Error\n";
    std::cout << std::string(85, '-') << "\n";
    
    for (double r : radii) {
        // Create data
        auto data = create_kerr_data_at_radius(r, theta, M, a);
        
        // For Kerr in Boyer-Lindquist coordinates, f(r) is derived from metric
        double Sigma = r*r + a*a*cos(theta)*cos(theta);
        double Delta = r*r - 2.0*M*r + a*a;
        double f = -a * (r*r + a*a - M*r) / Delta;
        
        // Compute angular momentum using axisymmetric formula
        double J_computed = (1.0/3.0) * data.alpha * r * r * f;
        
        // Error
        double error = J_computed - J_expected;
        double rel_error = std::abs(error) / std::abs(J_expected);
        
        std::cout << std::setw(10) << std::setprecision(1) << std::fixed << r
                  << std::setw(15) << std::setprecision(6) << std::scientific << data.alpha
                  << std::setw(15) << f
                  << std::setw(15) << J_computed
                  << std::setw(15) << error
                  << std::setw(15) << rel_error << "\n";
    }
    
    // Demonstrate vorticity bound
    std::cout << "\n\nVorticity Bound Demonstration:\n";
    std::cout << "Theorem 3.1: ||V ∧ dV||_g ≤ ||α||_g ||dβ||_g²\n";
    std::cout << std::string(50, '-') << "\n";
    
    double r = 100.0;
    auto data = create_kerr_data_at_radius(r, theta, M, a);
    
    auto [lhs, rhs] = compute_vorticity_bound(data);
    
    std::cout << "At r = " << r << ":\n";
    std::cout << "  ||V ∧ dV||_g = " << lhs << "\n";
    std::cout << "  ||α||_g ||dβ||_g² = " << rhs << "\n";
    std::cout << "  Bound satisfied: " << (lhs <= rhs ? "YES" : "NO") << "\n";
    std::cout << "  Ratio lhs/rhs = " << lhs/rhs << "\n";
    
    // Zero angular momentum condition check
    std::cout << "\nZero Angular Momentum Condition (Corollary 3.3):\n";
    std::cout << "If α ≡ 0 or dβ ≡ 0, then J = 0\n";
    
    std::vector<MetricData> test_data = {data};
    auto [alpha_zero, dbeta_zero] = check_zero_angular_momentum(test_data);
    
    std::cout << "  α ≡ 0? " << (alpha_zero ? "YES" : "NO") << "\n";
    std::cout << "  dβ ≡ 0? " << (dbeta_zero ? "YES" : "NO") << "\n";
    std::cout << "  ⇒ J = 0? " << ((alpha_zero || dbeta_zero) ? "YES" : "NO") << "\n";
    
    std::cout << "\n=============================================\n";
    std::cout << "Example completed successfully!\n";
    std::cout << "=============================================\n";
    
    return 0;
}