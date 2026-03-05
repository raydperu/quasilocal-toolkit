#include <iostream>
#include <vector>
#include <cmath>

// Forward declarations only - no headers with default arguments
namespace qlt {
    std::vector<double> compute_axisymmetric_j(
        const std::vector<double>& alpha,
        const std::vector<double>& f,
        const std::vector<double>& r);
        
    double compute_axisymmetric_surgical_flux(
        double alpha1, double r1, double f1,
        double alpha2, double r2, double f2);
}

int main() {
    std::cout << "=== QuasiLocal Toolkit Example ===\n\n";
    
    // Create test data
    std::vector<double> r = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> alpha = {0.9, 0.95, 0.98, 0.99, 1.0};
    std::vector<double> f = {0.1, 0.2, 0.3, 0.4, 0.5};
    
    std::cout << "Input data:\n";
    for (size_t i = 0; i < r.size(); i++) {
        std::cout << "  r = " << r[i] << ", α = " << alpha[i] << ", f = " << f[i] << "\n";
    }
    
    // Compute J(r)
    std::vector<double> J = qlt::compute_axisymmetric_j(alpha, f, r);
    
    std::cout << "\nResults J(r):\n";
    for (size_t i = 0; i < J.size(); i++) {
        std::cout << "  J(" << r[i] << ") = " << J[i] << "\n";
    }
    
    // Compute surgical flux
    double J_flux = qlt::compute_axisymmetric_surgical_flux(
        alpha[1], r[1], f[1],
        alpha[3], r[3], f[3]
    );
    std::cout << "\nSurgical flux between r=2 and r=4: " << J_flux << "\n";
    
    std::cout << "\n=== Example Complete ===\n";
    return 0;
}ECHO is off.
    double compute_axisymmetric_surgical_flux( 
        double alpha1, double r1, double f1, 
        double alpha2, double r2, double f2); 
} 
 
int main() { 
    std::cout << "=== QuasiLocal Toolkit Example ===\n\n"; 
    std::vector<double> r = {1.0, 2.0, 3.0, 4.0, 5.0}; 
    std::vector<double> alpha = {0.9, 0.95, 0.98, 0.99, 1.0}; 
    std::vector<double> f = {0.1, 0.2, 0.3, 0.4, 0.5}; 
    std::cout << "Input data:\n"; 
    for (size_t i = 0; i < r.size(); i++) { 
        std::cout << "  r = " << r[i] << ", ╬▒ = " << alpha[i] << ", f = " << f[i] << "\n"; 
    } 
    auto J = qlt::compute_axisymmetric_j(alpha, f, r); 
    std::cout << "\nResults J(r):\n"; 
    for (size_t i = 0; i < J.size(); i++) { 
        std::cout << "  J(" << r[i] << ") = " << J[i] << "\n"; 
    } 
    double flux = qlt::compute_axisymmetric_surgical_flux( 
        alpha[1], r[1], f[1], alpha[3], r[3], f[3]); 
    std::cout << "\nSurgical flux between r=2 and r=4: " << flux << "\n"; 
    std::cout << "\n=== Example Complete ===\n"; 
    return 0; 
} 
ECHO is off.
    double compute_axisymmetric_surgical_flux( 
        double alpha1, double r1, double f1, 
        double alpha2, double r2, double f2); 
} 
 
int main() { 
    std::cout << "=== QuasiLocal Toolkit Example ===\n\n"; 
    std::vector<double> r = {1.0, 2.0, 3.0, 4.0, 5.0}; 
    std::vector<double> alpha = {0.9, 0.95, 0.98, 0.99, 1.0}; 
    std::vector<double> f = {0.1, 0.2, 0.3, 0.4, 0.5}; 
    std::cout << "Input data:\n"; 
    for (size_t i = 0; i < r.size(); i++) { 
        std::cout << "  r = " << r[i] << ", ╬▒ = " << alpha[i] << ", f = " << f[i] << "\n"; 
    } 
    auto J = qlt::compute_axisymmetric_j(alpha, f, r); 
    std::cout << "\nResults J(r):\n"; 
    for (size_t i = 0; i < J.size(); i++) { 
        std::cout << "  J(" << r[i] << ") = " << J[i] << "\n"; 
    } 
    double flux = qlt::compute_axisymmetric_surgical_flux( 
        alpha[1], r[1], f[1], alpha[3], r[3], f[3]); 
    std::cout << "\nSurgical flux between r=2 and r=4: " << flux << "\n"; 
    std::cout << "\n=== Example Complete ===\n"; 
    return 0; 
} 
 
    double compute_axisymmetric_surgical_flux( 
        double alpha1, double r1, double f1, 
        double alpha2, double r2, double f2); 
} 
 
int main() { 
    std::cout << "=== QuasiLocal Toolkit Example ===\n\n"; 
 
    // Test data 
    std::vector<double> r = {1.0, 2.0, 3.0, 4.0, 5.0}; 
    std::vector<double> alpha = {0.9, 0.95, 0.98, 0.99, 1.0}; 
    std::vector<double> f = {0.1, 0.2, 0.3, 0.4, 0.5}; 
 
    std::cout << "Input data:\n"; 
    for (size_t i = 0; i < r.size(); i++) { 
        std::cout << "  r = " << r[i] 
                  << ", alpha = " << alpha[i] 
                  << ", f = " << f[i] << "\n"; 
    } 
 
    // Compute J(r) 
    std::vector<double> J = qlt::compute_axisymmetric_j(alpha, f, r); 
 
    std::cout << "\nResults J(r):\n"; 
    for (size_t i = 0; i < J.size(); i++) { 
        std::cout << "  J(" << r[i] << ") = " << J[i] << "\n"; 
    } 
 
    // Compute surgical flux between r=2 and r=4 
    double flux = qlt::compute_axisymmetric_surgical_flux( 
        alpha[1], r[1], f[1], 
        alpha[3], r[3], f[3]); 
    std::cout << "\nSurgical flux between r=2 and r=4: " << flux << "\n"; 
 
    std::cout << "\n=== Example Complete ===\n"; 
    return 0; 
} 
