// coordinate_systems.cpp 
#include "coordinate_systems.h" 
 
namespace qlt { 
namespace grid { 
 
    double r = std::sqrt(x*x + y*y + z*z); 
    double theta = std::acos(z / r); 
    double phi = std::atan2(y, x); 
    return {r, theta, phi}; 
} 
 
    double x = r * std::sin(theta) * std::cos(phi); 
    double y = r * std::sin(theta) * std::sin(phi); 
    double z = r * std::cos(theta); 
    return {x, y, z}; 
} 
 
    double rho = std::sqrt(x*x + y*y); 
    double phi = std::atan2(y, x); 
    return {rho, phi, z}; 
} 
 
    double x = rho * std::cos(phi); 
    double y = rho * std::sin(phi); 
    return {x, y, z}; 
} 
 
} // namespace grid 
} // namespace qlt 
