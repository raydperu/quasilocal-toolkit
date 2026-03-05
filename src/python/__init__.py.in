"""
QuasiLocal Toolkit (QLT) - Python interface

Mathematically accurate implementation of quasilocal conserved quantities
from "Geometric Vorticity and Quasilocal Angular Momentum Conservation"
by Rayan D. Peru.

This module provides Python bindings to the C++ quasilocal toolkit.
"""

__version__ = "@PROJECT_VERSION@"
__author__ = "Rayan D. Peru"
__license__ = "MIT"

# Import the compiled module
from .qlt import *

# Re-export main classes and functions
__all__ = [
    # Core classes
    'MetricData',
    
    # Core functions
    'compute_axisymmetric_j',
    'compute_clebsch_komar_j',
    'compute_quasilocal_mass',
    'compute_surgical_flux',
    'compute_vorticity_bound',
    'check_zero_angular_momentum',
    
    # Grid functions
    'MetricReader',
    'CoordinateSystems',
    'KillingVectorGenerator',
    'SurfaceIntegral',
    'VolumeIntegral',
    
    # Constants
    'M_PI',
]

# Add Pythonic convenience functions
def angular_momentum_sphere(radius, alpha, f):
    """
    Compute angular momentum on a sphere using axisymmetric formula.
    
    Parameters
    ----------
    radius : float or array-like
        Radius of the sphere
    alpha : float or array-like
        Clebsch potential α
    f : float or array-like
        Phase function f(r) = β - φ
    
    Returns
    -------
    J : float or array-like
        Angular momentum J(r) = (1/3) α r² f
    """
    import numpy as np
    
    radius = np.asarray(radius)
    alpha = np.asarray(alpha)
    f = np.asarray(f)
    
    return (1.0/3.0) * alpha * radius * radius * f

def validate_kerr_solution(M, a, r_values):
    """
    Validate against Kerr solution.
    
    Parameters
    ----------
    M : float
        Black hole mass
    a : float
        Spin parameter (a = J/M)
    r_values : array-like
        Radial coordinates to evaluate at
    
    Returns
    -------
    results : dict
        Dictionary with validation results
    """
    import numpy as np
    
    r = np.asarray(r_values)
    J_expected = M * a
    
    # Kerr α(r) in Boyer-Lindquist coordinates (simplified)
    # Actual implementation would use full Kerr metric
    alpha_kerr = -2.0 * M * a * r / (r * r + a * a)
    
    # Kerr f(r) satisfying Clebsch constraint
    f_kerr = -a * (r * r + a * a - M * r) / (r * r - 2.0 * M * r + a * a)
    
    J_computed = angular_momentum_sphere(r, alpha_kerr, f_kerr)
    
    # Relative error
    error = np.abs(J_computed - J_expected) / np.abs(J_expected)
    
    return {
        'r': r,
        'J_computed': J_computed,
        'J_expected': J_expected,
        'error': error,
        'max_error': np.max(error),
        'mean_error': np.mean(error)
    }

def create_test_sphere(num_theta=32, num_phi=64, radius=10.0):
    """
    Create a test spherical surface for demonstrations.
    
    Parameters
    ----------
    num_theta : int
        Number of θ points
    num_phi : int
        Number of φ points
    radius : float
        Sphere radius
    
    Returns
    -------
    sphere_data : list of MetricData
        List of MetricData objects on the sphere
    """
    import numpy as np
    
    sphere_data = []
    
    for i in range(num_theta):
        theta = np.pi * (i + 0.5) / num_theta  # Avoid poles
        
        for j in range(num_phi):
            phi = 2.0 * np.pi * j / num_phi
            
            # Create MetricData (simplified)
            md = MetricData()
            md.r = radius
            md.theta = theta
            md.phi_coord = phi
            
            # Flat metric in spherical coordinates
            md.h_ij = np.array([
                [1.0, 0.0, 0.0],
                [0.0, radius*radius, 0.0],
                [0.0, 0.0, radius*radius*np.sin(theta)*np.sin(theta)]
            ], dtype=np.float64)
            
            md.sqrt_det_h = np.sqrt(np.linalg.det(md.h_ij))
            
            # Rotational Killing vector about z-axis
            x = radius * np.sin(theta) * np.cos(phi)
            y = radius * np.sin(theta) * np.sin(phi)
            md.xi = np.array([-y, x, 0.0], dtype=np.float64)
            
            sphere_data.append(md)
    
    return sphere_data

# Add documentation
__doc__ = __doc__ + """

Examples
--------
>>> import qlt
>>> import numpy as np
>>>
>>> # Compute angular momentum using axisymmetric formula
>>> r = np.linspace(1.0, 100.0, 100)
>>> alpha = -0.1 * r  # Example α(r)
>>> f = 0.5 / r       # Example f(r)
>>> J = qlt.angular_momentum_sphere(r, alpha, f)
>>>
>>> # Validate against Kerr solution
>>> results = qlt.validate_kerr_solution(M=1.0, a=0.5, r_values=r)
>>> print(f"Maximum error: {results['max_error']:.2e}")
>>>
>>> # Create a test sphere
>>> sphere = qlt.create_test_sphere(num_theta=16, num_phi=32, radius=10.0)
>>> print(f"Created sphere with {len(sphere)} points")
"""

# Version info
def version():
    """Return package version."""
    return __version__

def citation():
    """Return citation information."""
    return """Rayan D. Peru, "Geometric Vorticity and Quasilocal Angular Momentum Conservation", 2025."""