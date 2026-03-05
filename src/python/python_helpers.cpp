#include "bindings.h"
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace qlt {
namespace python {

// Python-specific helper functions

/**
 * @brief Convert Python list to std::vector<MetricData>
 */
std::vector<MetricData> py_list_to_metric_data(const py::list &py_list) {
    std::vector<MetricData> result;
    result.reserve(py_list.size());
    
    for (const auto &item : py_list) {
        if (py::isinstance<MetricData>(item)) {
            result.push_back(item.cast<MetricData>());
        } else if (py::isinstance<py::dict>(item)) {
            // Convert from dict
            py::dict d = item.cast<py::dict>();
            MetricData md;
            
            if (d.contains("h_ij")) md.h_ij = d["h_ij"].cast<Matrix3d>();
            if (d.contains("K_ij")) md.K_ij = d["K_ij"].cast<Matrix3d>();
            if (d.contains("xi")) md.xi = d["xi"].cast<Vector3d>();
            if (d.contains("alpha")) md.alpha = d["alpha"].cast<double>();
            if (d.contains("phi")) md.phi = d["phi"].cast<double>();
            if (d.contains("beta")) md.beta = d["beta"].cast<double>();
            if (d.contains("r")) md.r = d["r"].cast<double>();
            if (d.contains("theta")) md.theta = d["theta"].cast<double>();
            if (d.contains("phi_coord")) md.phi_coord = d["phi_coord"].cast<double>();
            
            result.push_back(md);
        }
    }
    
    return result;
}

/**
 * @brief Convert std::vector<MetricData> to Python list
 */
py::list metric_data_to_py_list(const std::vector<MetricData> &data) {
    py::list result;
    
    for (const auto &md : data) {
        result.append(md);
    }
    
    return result;
}

/**
 * @brief Create test data for Kerr metric
 */
py::dict create_kerr_test_data(double M, double a, 
                               const py::array_t<double> &r_values,
                               double theta = M_PI/2.0) {
    auto r_buf = r_values.request();
    size_t n = r_buf.size;
    
    std::vector<double> r_vec(n);
    std::vector<double> alpha_vec(n);
    std::vector<double> f_vec(n);
    std::vector<double> J_vec(n);
    
    // Copy r values
    std::memcpy(r_vec.data(), r_buf.ptr, n * sizeof(double));
    
    // Kerr metric in Boyer-Lindquist coordinates
    for (size_t i = 0; i < n; i++) {
        double r = r_vec[i];
        double Sigma = r*r + a*a*cos(theta)*cos(theta);
        double Delta = r*r - 2.0*M*r + a*a;
        
        // Kerr α(r)
        alpha_vec[i] = -2.0 * M * a * r * sin(theta)*sin(theta) / Sigma;
        
        // f(r) satisfying Clebsch constraint
        f_vec[i] = -a * (r*r + a*a - M*r) / Delta;
        
        // Angular momentum
        J_vec[i] = (1.0/3.0) * alpha_vec[i] * r * r * f_vec[i];
    }
    
    py::dict result;
    result["r"] = r_values;
    result["alpha"] = py::array_t<double>(n, alpha_vec.data());
    result["f"] = py::array_t<double>(n, f_vec.data());
    result["J"] = py::array_t<double>(n, J_vec.data());
    result["J_expected"] = M * a;
    
    return result;
}

/**
 * @brief Create test data for Schwarzschild metric
 */
py::dict create_schwarzschild_test_data(double M,
                                        const py::array_t<double> &r_values) {
    auto r_buf = r_values.request();
    size_t n = r_buf.size;
    
    std::vector<double> r_vec(n);
    std::vector<double> mass_vec(n);
    
    std::memcpy(r_vec.data(), r_buf.ptr, n * sizeof(double));
    
    for (size_t i = 0; i < n; i++) {
        double r = r_vec[i];
        // For Schwarzschild, mass aspect should be M at all r > 2M
        mass_vec[i] = M;
    }
    
    py::dict result;
    result["r"] = r_values;
    result["M"] = py::array_t<double>(n, mass_vec.data());
    result["M_expected"] = M;
    
    return result;
}

/**
 * @brief Register Python helper functions
 */
void bind_python_helpers(py::module &m) {
    m.def("create_kerr_test_data", &create_kerr_test_data,
        py::arg("M"), py::arg("a"), py::arg("r_values"),
        py::arg("theta") = M_PI/2.0,
        "Create test data for Kerr metric validation")
    
    .def("create_schwarzschild_test_data", &create_schwarzschild_test_data,
        py::arg("M"), py::arg("r_values"),
        "Create test data for Schwarzschild metric validation")
    
    .def("py_list_to_metric_data", &py_list_to_metric_data,
        py::arg("py_list"),
        "Convert Python list to vector of MetricData")
    
    .def("metric_data_to_py_list", &metric_data_to_py_list,
        py::arg("data"),
        "Convert vector of MetricData to Python list");
}

} // namespace python
} // namespace qlt