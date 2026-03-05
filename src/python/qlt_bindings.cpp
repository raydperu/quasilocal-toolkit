#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <vector>
#include <array>
#include <Eigen/Dense>

namespace py = pybind11;

// Forward declare the namespace and functions
namespace qlt {
    // Complete definition of MetricData
    struct MetricData {
        Eigen::Matrix3d h_ij;           // Spatial metric
        Eigen::Matrix3d K_ij;            // Extrinsic curvature
        Eigen::Vector3d xi;               // Killing vector field
        Eigen::Matrix3d dxi_dx;           // Derivative ∂_i ξ_j
        double sqrt_det_h;                 // √det(h_ij)
        double alpha;                      // Clebsch potential α
        double phi;                         // Clebsch potential Φ
        double beta;                        // Clebsch potential β
        Eigen::Vector3d dalpha_dx;          // ∇α
        Eigen::Vector3d dbeta_dx;           // ∇β
        Eigen::Vector3d dphi_dx;            // ∇Φ
        Eigen::Vector3d normal;             // Unit normal
        double r;                            // Radial coordinate
        double theta;                        // Polar angle
        double phi_coord;                    // Azimuthal angle
        
        // Default constructor
        MetricData() 
            : h_ij(Eigen::Matrix3d::Identity())
            , K_ij(Eigen::Matrix3d::Zero())
            , xi(Eigen::Vector3d::Zero())
            , dxi_dx(Eigen::Matrix3d::Zero())
            , sqrt_det_h(1.0)
            , alpha(0.0)
            , phi(0.0)
            , beta(0.0)
            , dalpha_dx(Eigen::Vector3d::Zero())
            , dbeta_dx(Eigen::Vector3d::Zero())
            , dphi_dx(Eigen::Vector3d::Zero())
            , normal(Eigen::Vector3d::Zero())
            , r(0.0)
            , theta(0.0)
            , phi_coord(0.0) {}
    };
    
    // Function declarations
    std::vector<double> compute_axisymmetric_j(
        const std::vector<double>& alpha,
        const std::vector<double>& f,
        const std::vector<double>& r);
    
    double compute_axisymmetric_surgical_flux(
        double alpha1, double r1, double f1,
        double alpha2, double r2, double f2);
    
    std::vector<double> extract_alpha_from_surface(
        const std::vector<MetricData>& surface_data);
    
    std::vector<double> extract_f_from_surface(
        const std::vector<MetricData>& surface_data);
    
    std::vector<double> extract_r_from_surface(
        const std::vector<MetricData>& surface_data);
    
    double compute_clebsch_komar_j(
        const std::vector<MetricData>& surface_data,
        const std::vector<MetricData>& volume_data,
        double kappa);
    
    double compute_quasilocal_mass(
        const std::vector<MetricData>& surface_data,
        bool use_reference);
    
    double compute_surgical_flux(
        const std::vector<MetricData>& inner_data,
        const std::vector<MetricData>& outer_data,
        bool axisymmetric);
    
    std::pair<bool, bool> check_zero_angular_momentum(
        const std::vector<MetricData>& data,
        double tolerance);
    
    std::pair<double, double> compute_vorticity_bound(
        const MetricData& data);
}

using namespace qlt;

// Forward declarations of binding functions
void init_axisymmetric_j(py::module_&);
void init_clebsch_komar(py::module_&);
void init_quasilocal_mass(py::module_&);
void init_surgical_flux(py::module_&);
void init_vorticity_bounds(py::module_&);

PYBIND11_MODULE(qlt, m) {
    m.doc() = "QuasiLocal Toolkit Python bindings";
    
    // Bind MetricData first so it's available for other bindings
    py::class_<MetricData>(m, "MetricData")
        .def(py::init<>())
        .def_readwrite("r", &MetricData::r)
        .def_readwrite("theta", &MetricData::theta)
        .def_readwrite("phi_coord", &MetricData::phi_coord)
        .def_readwrite("alpha", &MetricData::alpha)
        .def_readwrite("phi", &MetricData::phi)
        .def_readwrite("beta", &MetricData::beta)
        .def_readwrite("xi", &MetricData::xi)
        .def_readwrite("normal", &MetricData::normal)
        .def_readwrite("sqrt_det_h", &MetricData::sqrt_det_h)
        .def_property("h_ij", 
            [](const MetricData& d) { return d.h_ij; },
            [](MetricData& d, const Eigen::Matrix3d& v) { d.h_ij = v; })
        .def_property("K_ij",
            [](const MetricData& d) { return d.K_ij; },
            [](MetricData& d, const Eigen::Matrix3d& v) { d.K_ij = v; })
        .def_property("dxi_dx",
            [](const MetricData& d) { return d.dxi_dx; },
            [](MetricData& d, const Eigen::Matrix3d& v) { d.dxi_dx = v; })
        .def_property("dalpha_dx",
            [](const MetricData& d) { return d.dalpha_dx; },
            [](MetricData& d, const Eigen::Vector3d& v) { d.dalpha_dx = v; })
        .def_property("dbeta_dx",
            [](const MetricData& d) { return d.dbeta_dx; },
            [](MetricData& d, const Eigen::Vector3d& v) { d.dbeta_dx = v; })
        .def_property("dphi_dx",
            [](const MetricData& d) { return d.dphi_dx; },
            [](MetricData& d, const Eigen::Vector3d& v) { d.dphi_dx = v; })
        .def("__repr__", [](const MetricData& d) {
            return "<qlt.MetricData r=" + std::to_string(d.r) + 
                   " theta=" + std::to_string(d.theta) + ">";
        });
    
    // Initialize all submodules
    init_axisymmetric_j(m);
    init_clebsch_komar(m);
    init_quasilocal_mass(m);
    init_surgical_flux(m);
    init_vorticity_bounds(m);
    
    // Version information
    m.attr("__version__") = "0.1.0";
}

// Axisymmetric J bindings
void init_axisymmetric_j(py::module_& m) {
    py::module_ sub = m.def_submodule("axisymmetric_j", "Axisymmetric J computations");
    
    sub.def("compute", &compute_axisymmetric_j, 
            "Compute axisymmetric angular momentum J(r) = (1/3) α(r) r² f(r)",
            py::arg("alpha"), py::arg("f"), py::arg("r"));
    
    sub.def("surgical_flux", &compute_axisymmetric_surgical_flux,
            "Compute axisymmetric surgical flux ΔJ = (1/3)[α(r₂)r₂²f(r₂) - α(r₁)r₁²f(r₁)]",
            py::arg("alpha1"), py::arg("r1"), py::arg("f1"),
            py::arg("alpha2"), py::arg("r2"), py::arg("f2"));
    
    sub.def("extract_alpha", &extract_alpha_from_surface,
            "Extract α(r) from axisymmetric data",
            py::arg("surface_data"));
    
    sub.def("extract_f", &extract_f_from_surface,
            "Extract f(r) = β - φ from axisymmetric data",
            py::arg("surface_data"));
    
    sub.def("extract_r", &extract_r_from_surface,
            "Extract radial coordinates from surface data",
            py::arg("surface_data"));
}

// Clebsch-Komar bindings
void init_clebsch_komar(py::module_& m) {
    py::module_ sub = m.def_submodule("clebsch_komar", "Clebsch-Komar quantities");
    
    sub.def("compute", &compute_clebsch_komar_j,
            "Compute Clebsch-Komar angular momentum J_S[ξ]",
            py::arg("surface_data"), py::arg("volume_data"), py::arg("kappa"));
}

// Quasilocal mass bindings
void init_quasilocal_mass(py::module_& m) {
    py::module_ sub = m.def_submodule("quasilocal_mass", "Quasilocal mass computations");
    
    sub.def("compute", &compute_quasilocal_mass,
            "Compute quasilocal mass via modified Brown-York formalism",
            py::arg("surface_data"), py::arg("use_reference"));
}

// Surgical flux bindings
void init_surgical_flux(py::module_& m) {
    py::module_ sub = m.def_submodule("surgical_flux", "Surgical flux computations");
    
    sub.def("compute", &compute_surgical_flux,
            "Compute surgical flux across excised region",
            py::arg("inner_data"), py::arg("outer_data"), py::arg("axisymmetric"));
}

// Vorticity bounds bindings
void init_vorticity_bounds(py::module_& m) {
    py::module_ sub = m.def_submodule("vorticity_bounds", "Vorticity bounds computations");
    
    sub.def("check_zero_angular_momentum", &check_zero_angular_momentum,
            "Check if zero angular momentum condition is satisfied",
            py::arg("data"), py::arg("tolerance"));
    
    sub.def("compute_bound", &compute_vorticity_bound,
            "Compute vorticity bound ||V ∧ dV||_g ≤ ||α||_g ||dβ||_g²",
            py::arg("data"));
}