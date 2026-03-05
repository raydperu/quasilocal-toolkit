#ifndef QLT_COORDINATE_SYSTEMS_H 
#define QLT_COORDINATE_SENSITIVITY_H 
 
#include "../core/quasilocal_quantities.h" 
#include <vector> 
#include <array> 
#include <cmath> 
#include <string> 
#include <stdexcept> 
 
namespace qlt { 
namespace grid { 
 
class CoordinateSystems { 
public: 
    static std::array<double, 3> cartesian_to_spherical(double x, double y, double z); 
    static std::array<double, 3> spherical_to_cartesian(double r, double theta, double phi); 
}; 
 
} // namespace grid 
} // namespace qlt 
 
#endif // QLT_COORDINATE_SYSTEMS_H 
