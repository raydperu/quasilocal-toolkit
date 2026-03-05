#ifndef QLT_PYTHON_BINDINGS_H
#define QLT_PYTHON_BINDINGS_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "qlt/quasilocal_quantities.h"
#include "grid/metric_reader.h"
#include "grid/coordinate_systems.h"
#include "grid/killing_vector.h"
#include "grid/surface_integral.h"
#include "grid/volume_integral.h"

namespace py = pybind11;

// Forward declarations
void bind_core_module(py::module &m);
void bind_grid_module(py::module &m);
void bind_metric_data(py::module &m);
void bind_enums(py::module &m);

#endif // QLT_PYTHON_BINDINGS_H