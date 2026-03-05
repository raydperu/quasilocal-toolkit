# QuasiLocal Toolkit Documentation

## Overview

The **QuasiLocal Toolkit (QLT)** is a mathematically accurate implementation
of quasilocal conserved quantities in general relativity, based on the
framework developed in:

> **"Geometric Vorticity and Quasilocal Angular Momentum Conservation"**  
> by Rayan D. Peru

This toolkit provides computational tools for tracking angular momentum,
mass, and other conserved quantities in numerical relativity simulations
without requiring horizon finding algorithms.

## Key Features

### Mathematical Accuracy
- **Exact implementation** of formulas from the paper
- **Rigorous validation** against analytic solutions (Kerr, Schwarzschild)
- **Mathematical proofs** translated to computational algorithms

### Computational Efficiency
- **No horizon finding required** - uses Clebsch field values directly
- **Axisymmetric formula**: \(J(r) = \frac{1}{3}\alpha(r)r^2f(r)\)
- **Surgical flux computation** for excision regions

### Numerical Relativity Integration
- **Data I/O** for SpEC, Cactus HDF5/ASCII formats
- **Coordinate transformations** between common systems
- **Killing vector generation** for common spacetimes

### Production Quality
- **Comprehensive test suite** with Catch2
- **Python bindings** with PyBind11
- **Modern CMake build system**
- **API documentation** with Doxygen

## Mathematical Framework

### Core Formulas

#### 1. Axisymmetric Angular Momentum (Theorem 4.3)
\[
J(r) = \frac{1}{3}\alpha(r)r^2f(r)
\]
where:
- \(\alpha(r)\) = Clebsch amplitude field
- \(f(r) = \beta - \phi\) = phase function
- \(r\) = radial coordinate

#### 2. Clebsch-Komar Identity (Eq. 1.3)
\[
J_S[\xi] = \frac{1}{16\pi}\int_S \nabla^\mu \xi^\nu dS_{\mu\nu} 
         + \frac{\kappa}{16\pi}\int_{\Sigma_S} [\alpha d\beta \wedge d\xi^\flat 
         + \Phi d(d\beta \wedge \xi^\flat)]
\]

#### 3. Surgical Flux Theorem (Theorem 4.4)
\[
\Delta J = J(\lambda^+) - J(\lambda^-) = -\frac{1}{8\pi}\int_{\partial V} \alpha d\beta \wedge \xi^\flat
\]

#### 4. Vorticity Bounds (Theorem 3.1/4.1)
\[
\|V \wedge dV\|_g \leq \|\alpha\|_g \|d\beta\|_g^2
\]

#### 5. Zero Angular Momentum Condition (Corollary 3.3/4.3)
\[
\alpha \equiv 0 \quad \text{or} \quad d\beta \equiv 0 \quad \Rightarrow \quad J = 0
\]

### Physical Interpretation

The Clebsch representation \(V_\mu = \nabla_\mu\Phi + \alpha\nabla_\mu\beta\):
- **\(\Phi\)**: Velocity potential (irrotational part)
- **\(\alpha\)**: Amplitude of rotational structure
- **\(\beta\)**: Phase of rotational structure
- **\(V_\mu\)**: Vorticity 1-form (geometric analogue of fluid vorticity)

## Getting Started

### Installation

```bash
# Clone repository
git clone https://github.com/username/quasilocal-toolkit.git
cd quasilocal-toolkit

# Build
mkdir build && cd build
cmake -DBUILD_PYTHON_BINDINGS=ON -DBUILD_TESTS=ON ..
make -j4

# Test
ctest

# Install
make install
# or for Python
pip install .