# QuasiLocal Toolkit Examples

This directory contains example programs demonstrating the use of the
QuasiLocal Toolkit for computing quasilocal conserved quantities.

## C++ Examples

### 1. Kerr Validation (`kerr_validation.cpp`)
Validates the quasilocal angular momentum computation against the Kerr metric.

**Key features:**
- Creates Kerr metric data at various radii
- Computes angular momentum using axisymmetric formula
- Compares with expected value J = M*a
- Demonstrates vorticity bounds
- Checks zero angular momentum condition

**Build and run:**
```bash
cd build
cmake -DBUILD_EXAMPLES=ON ..
make example_kerr
./examples/example_kerr