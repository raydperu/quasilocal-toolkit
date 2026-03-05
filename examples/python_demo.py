#!/usr/bin/env python3
"""
Python demonstration of quasilocal toolkit.
"""

import numpy as np
import matplotlib.pyplot as plt

def demo_axisymmetric_formula():
    """Demonstrate axisymmetric angular momentum formula."""
    print("=" * 60)
    print("Axisymmetric Angular Momentum Formula")
    print("=" * 60)
    
    # Create test data
    r = np.linspace(1.0, 100.0, 100)
    alpha = -0.1 * r  # Linear α(r)
    f = 0.5 / r       # f(r) ~ 1/r
    
    # Compute angular momentum
    J = (1.0/3.0) * alpha * r * r * f
    
    print(f"Radii: {len(r)} points from {r[0]} to {r[-1]}")
    print(f"α(r) range: [{alpha.min():.3f}, {alpha.max():.3f}]")
    print(f"f(r) range: [{f.min():.3f}, {f.max():.3f}]")
    print(f"J(r) range: [{J.min():.3f}, {J.max():.3f}]")
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    
    axes[0, 0].plot(r, alpha)
    axes[0, 0].set_xlabel('r')
    axes[0, 0].set_ylabel('α(r)')
    axes[0, 0].set_title('Clebsch potential α(r)')
    axes[0, 0].grid(True)
    
    axes[0, 1].plot(r, f)
    axes[0, 1].set_xlabel('r')
    axes[0, 1].set_ylabel('f(r)')
    axes[0, 1].set_title('Phase function f(r) = β - φ')
    axes[0, 1].grid(True)
    
    axes[1, 0].plot(r, J)
    axes[1, 0].set_xlabel('r')
    axes[1, 0].set_ylabel('J(r)')
    axes[1, 0].set_title('Angular momentum J(r) = (1/3) α r² f')
    axes[1, 0].grid(True)
    
    # Show product α r² f
    product = alpha * r * r * f
    axes[1, 1].plot(r, product)
    axes[1, 1].set_xlabel('r')
    axes[1, 1].set_ylabel('α r² f')
    axes[1, 1].set_title('Product α r² f')
    axes[1, 1].grid(True)
    
    plt.tight_layout()
    plt.savefig('axisymmetric_demo.png', dpi=150)
    print(f"\nPlot saved as 'axisymmetric_demo.png'")
    
    return J

def demo_kerr_validation():
    """Validate against Kerr solution."""
    print("\n" + "=" * 60)
    print("Kerr Metric Validation")
    print("=" * 60)
    
    # Kerr parameters
    M = 1.0      # Mass
    a = 0.5      # Spin parameter
    J_expected = M * a
    
    # Test radii
    r = np.logspace(0.5, 3, 50)  # Logarithmic spacing
    
    # Kerr metric in Boyer-Lindquist coordinates
    theta = np.pi/2.0  # Equatorial plane
    Sigma = r*r + a*a*np.cos(theta)*np.cos(theta)
    Delta = r*r - 2.0*M*r + a*a
    
    # Kerr Clebsch potentials
    alpha_kerr = -2.0 * M * a * r * np.sin(theta)*np.sin(theta) / Sigma
    f_kerr = -a * (r*r + a*a - M*r) / Delta
    
    # Compute angular momentum
    J_computed = (1.0/3.0) * alpha_kerr * r * r * f_kerr
    
    # Error
    error = np.abs(J_computed - J_expected)
    rel_error = error / np.abs(J_expected)
    
    print(f"Kerr parameters: M = {M}, a = {a}")
    print(f"Expected angular momentum: J = M*a = {J_expected}")
    print(f"Maximum relative error: {rel_error.max():.2e}")
    print(f"Mean relative error: {rel_error.mean():.2e}")
    
    # Plot convergence
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    axes[0].loglog(r, error)
    axes[0].set_xlabel('Radius r')
    axes[0].set_ylabel('Absolute error |J - M*a|')
    axes[0].set_title('Convergence to Kerr angular momentum')
    axes[0].grid(True, which='both')
    
    axes[1].loglog(r, rel_error)
    axes[1].set_xlabel('Radius r')
    axes[1].set_ylabel('Relative error |J - M*a|/|M*a|')
    axes[1].set_title('Relative error')
    axes[1].grid(True, which='both')
    
    # Add power law fits
    for i in [0, 1]:
        # Fit power law to last 10 points
        mask = r > 50
        if mask.sum() > 5:
            coeffs = np.polyfit(np.log(r[mask]), np.log(error[mask]), 1)
            slope = coeffs[0]
            axes[i].text(0.05, 0.95, f'Slope: {slope:.2f}', 
                        transform=axes[i].transAxes, verticalalignment='top')
    
    plt.tight_layout()
    plt.savefig('kerr_validation.png', dpi=150)
    print(f"\nPlot saved as 'kerr_validation.png'")
    
    return J_computed, error

def demo_surgical_flux():
    """Demonstrate surgical flux computation."""
    print("\n" + "=" * 60)
    print("Surgical Flux Demonstration")
    print("=" * 60)
    
    # Parameters for surgical excision
    r_inner = 1.0
    r_outer = 2.0
    alpha_inner = 2.0
    alpha_outer = 1.0
    f_inner = 0.5
    f_outer = 0.25
    
    print(f"Surgical region: r₁ = {r_inner}, r₂ = {r_outer}")
    print(f"Clebsch potentials: α₁ = {alpha_inner}, α₂ = {alpha_outer}")
    print(f"Phase functions: f₁ = {f_inner}, f₂ = {f_outer}")
    
    # Axisymmetric surgical flux formula
    term_inner = alpha_inner * r_inner * r_inner * f_inner
    term_outer = alpha_outer * r_outer * r_outer * f_outer
    delta_J = (1.0/3.0) * (term_outer - term_inner)
    
    print(f"\nAxisymmetric surgical flux (Corollary 4.6):")
    print(f"  ΔJ = (1/3)[α₂r₂²f₂ - α₁r₁²f₁]")
    print(f"     = (1/3)[{alpha_outer}×{r_outer}²×{f_outer} - {alpha_inner}×{r_inner}²×{f_inner}]")
    print(f"     = (1/3)[{term_outer} - {term_inner}]")
    print(f"     = {delta_J}")
    
    # Demonstrate with multiple surgeries
    print(f"\nMultiple surgical excisions:")
    
    radii = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    alphas = np.array([2.0, 1.5, 1.0, 0.5, 0.0])
    fs = 0.5 / radii
    
    cumulative_loss = 0.0
    for i in range(len(radii) - 1):
        r1, r2 = radii[i], radii[i+1]
        a1, a2 = alphas[i], alphas[i+1]
        f1, f2 = fs[i], fs[i+1]
        
        term1 = a1 * r1 * r1 * f1
        term2 = a2 * r2 * r2 * f2
        delta = (1.0/3.0) * (term2 - term1)
        cumulative_loss += delta
        
        print(f"  Surgery {i+1}: r={r1:.1f}→{r2:.1f}, ΔJ={delta:.3f}, cumulative={cumulative_loss:.3f}")
    
    print(f"\nTotal angular momentum loss: {cumulative_loss:.3f}")
    
    # Plot angular momentum profile
    J_profile = (1.0/3.0) * alphas * radii * radii * fs
    
    plt.figure(figsize=(10, 6))
    plt.plot(radii, J_profile, 'o-', linewidth=2, markersize=8)
    
    # Mark surgical boundaries
    for i in range(len(radii) - 1):
        plt.axvline(x=radii[i+1], color='red', linestyle='--', alpha=0.5)
        plt.text(radii[i+1], J_profile.max() * 0.9, f'Surgery {i+1}', 
                rotation=90, verticalalignment='top')
    
    plt.xlabel('Radius r')
    plt.ylabel('Angular momentum J(r)')
    plt.title('Angular momentum profile with surgical excisions')
    plt.grid(True)
    plt.savefig('surgical_flux_demo.png', dpi=150)
    print(f"\nPlot saved as 'surgical_flux_demo.png'")
    
    return delta_J, cumulative_loss

def main():
    """Run all demonstrations."""
    print("QuasiLocal Toolkit Python Demonstration")
    print("=" * 60)
    
    # Try to import the toolkit
    try:
        import qlt
        print(f"✓ Imported qlt version {qlt.version()}")
        print(f"  {qlt.citation()}")
    except ImportError:
        print("✗ Could not import qlt module")
        print("  Make sure the toolkit is built and installed")
        return
    
    # Run demonstrations
    demo_axisymmetric_formula()
    demo_kerr_validation()
    demo_surgical_flux()
    
    print("\n" + "=" * 60)
    print("All demonstrations completed!")
    print("Check the generated PNG files for plots.")
    print("=" * 60)

if __name__ == "__main__":
    main()