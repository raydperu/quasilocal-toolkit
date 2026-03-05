#ifndef QLT_PYTHON_HELPERS_H
#define QLT_PYTHON_HELPERS_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "qlt/quasilocal_quantities.h"

namespace py = pybind11;

namespace qlt {
namespace python {

// MetricData conversions
std::vector<MetricData> py_list_to_metric_data(const py::list &py_list);
py::list metric_data_to_py_list(const std::vector<MetricData> &data);

// Test data creation
py::dict create_kerr_test_data(double M, double a, 
                               const py::array_t<double> &r_values,
                               double theta = M_PI/2.0);

py::dict create_schwarzschild_test_data(double M,
                                        const py::array_t<double> &r_values);

// Register Python helper functions
void bind_python_helpers(py::module &m);

} // namespace python
} // namespace qlt

#endif // QLT_PYTHON_HELPERS_H