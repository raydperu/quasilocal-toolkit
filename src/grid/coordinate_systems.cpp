#include "coordinate_systems.h" 
#include <cmath> 
 
namespace qlt { 
namespace grid { 
 
std::array<double, 3> CoordinateSystems::cartesian_to_spherical(double x, double y, double z) { 
    double r = std::sqrt(x*x + y*y + z*z); 
    double theta = std::acos(z / r); 
    double phi = std::atan2(y, x); 
    return {r, theta, phi}; 
} 
 
std::array<double, 3> CoordinateSystems::spherical_to_cartesian(double r, double theta, double phi) { 
    double x = r * std::sin(theta) * std::cos(phi); 
    double y = r * std::sin(theta) * std::sin(phi); 
    double z = r * std::cos(theta); 
    return {x, y, z}; 
} 
 
} // namespace grid 
} // namespace qlt 
