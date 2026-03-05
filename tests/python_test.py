#!/usr/bin/env python3
"""
Python tests for quasilocal toolkit.
"""

import sys
import numpy as np

def test_import():
    """Test basic import."""
    try:
        import qlt
        print("✓ Import successful")
        print(f"  Version: {qlt.version()}")
        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False

def test_axisymmetric_formula():
    """Test axisymmetric angular momentum formula."""
    import qlt
    
    # Simple test
    r = np.array([1.0, 2.0, 3.0])
    alpha = np.array([1.0, 2.0, 3.0])
    f = np.array([1.0, 1.0, 1.0])
    
    J_expected = (1.0/3.0) * alpha * r * r * f
    J_computed = qlt.angular_momentum_sphere(r, alpha, f)
    
    if np.allclose(J_computed, J_expected):
        print("✓ Axisymmetric formula test passed")
        return True
    else:
        print("✗ Axisymmetric formula test failed")
        print(f"  Expected: {J_expected}")
        print(f"  Computed: {J_computed}")
        return False

def test_kerr_validation():
    """Test Kerr validation function."""
    import qlt
    
    M = 1.0
    a = 0.5
    r_values = np.linspace(10.0, 100.0, 10)
    
    results = qlt.validate_kerr_solution(M, a, r_values)
    
    if results['max_error'] < 0.1:  # 10% tolerance
        print(f"✓ Kerr validation passed (max error: {results['max_error']:.2e})")
        return True
    else:
        print(f"✗ Kerr validation failed (max error: {results['max_error']:.2e})")
        return False

def test_create_test_sphere():
    """Test sphere creation."""
    import qlt
    
    sphere = qlt.create_test_sphere(num_theta=8, num_phi=16, radius=5.0)
    
    if len(sphere) == 8 * 16:
        print(f"✓ Created sphere with {len(sphere)} points")
        
        # Check some properties
        first_point = sphere[0]
        if hasattr(first_point, 'r') and first_point.r == 5.0:
            print("  First point has correct radius")
            return True
        else:
            print("  First point properties incorrect")
            return False
    else:
        print(f"✗ Sphere creation failed: expected 128 points, got {len(sphere)}")
        return False

def main():
    """Run all tests."""
    print("Running Python tests for QuasiLocal Toolkit")
    print("=" * 50)
    
    tests = [
        test_import,
        test_axisymmetric_formula,
        test_kerr_validation,
        test_create_test_sphere,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"✗ {test.__name__} raised exception: {e}")
            failed += 1
    
    print("=" * 50)
    print(f"Results: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("All tests passed!")
        return 0
    else:
        print("Some tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())