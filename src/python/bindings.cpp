#include "bindings.h"
#include <iostream>
#include <sstream>

namespace qlt {
namespace python {

// ============================================================================
// MetricData Python binding
// ============================================================================

void bind_metric_data(py::module &m) {
    py::class_<MetricData>(m, "MetricData", "Metric and field data at a point")
        .def(py::init<>())
        
        // Metric data
        .def_readwrite("h_ij", &MetricData::h_ij, "Spatial metric (3x3)")
        .def_readwrite("K_ij", &MetricData::K_ij, "Extrinsic curvature (3x3)")
        .def_readwrite("sqrt_det_h", &MetricData::sqrt_det_h, "√det(h_ij)")
        
        // Killing vector
        .def_readwrite("xi", &MetricData::xi, "Killing vector field")
        .def_readwrite("dxi_dx", &MetricData::dxi_dx, "Derivative of Killing vector")
        
        // Clebsch potentials
        .def_readwrite("alpha", &MetricData::alpha, "Clebsch potential α")
        .def_readwrite("phi", &MetricData::phi, "Clebsch potential Φ")
        .def_readwrite("beta", &MetricData::beta, "Clebsch potential β")
        
        // Derivatives
        .def_readwrite("dalpha_dx", &MetricData::dalpha_dx, "∇α")
        .def_readwrite("dbeta_dx", &MetricData::dbeta_dx, "∇β")
        .def_readwrite("dphi_dx", &MetricData::dphi_dx, "∇Φ")
        
        // Position
        .def_readwrite("r", &MetricData::r, "Radial coordinate")
        .def_readwrite("theta", &MetricData::theta, "Polar angle")
        .def_readwrite("phi_coord", &MetricData::phi_coord, "Azimuthal angle")
        
        // Normal vector
        .def_readwrite("normal", &MetricData::normal, "Unit normal vector")
        
        // String representation
        .def("__repr__", [](const MetricData &d) {
            std::ostringstream oss;
            oss << "MetricData(r=" << d.r 
                << ", θ=" << d.theta 
                << ", φ=" << d.phi_coord
                << ", α=" << d.alpha
                << ", |ξ|=" << d.xi.norm() << ")";
            return oss.str();
        })
        
        // Convert to/from Python dict
        .def("to_dict", [](const MetricData &d) {
            py::dict dict;
            dict["h_ij"] = d.h_ij;
            dict["K_ij"] = d.K_ij;
            dict["xi"] = d.xi;
            dict["alpha"] = d.alpha;
            dict["phi"] = d.phi;
            dict["beta"] = d.beta;
            dict["r"] = d.r;
            dict["theta"] = d.theta;
            dict["phi_coord"] = d.phi_coord;
            return dict;
        }, "Convert to Python dictionary")
        
        .def_static("from_dict", [](const py::dict &dict) {
            MetricData d;
            if (dict.contains("h_ij")) d.h_ij = dict["h_ij"].cast<Matrix3d>();
            if (dict.contains("K_ij")) d.K_ij = dict["K_ij"].cast<Matrix3d>();
            if (dict.contains("xi")) d.xi = dict["xi"].cast<Vector3d>();
            if (dict.contains("alpha")) d.alpha = dict["alpha"].cast<double>();
            if (dict.contains("phi")) d.phi = dict["phi"].cast<double>();
            if (dict.contains("beta")) d.beta = dict["beta"].cast<double>();
            if (dict.contains("r")) d.r = dict["r"].cast<double>();
            if (dict.contains("theta")) d.theta = dict["theta"].cast<double>();
            if (dict.contains("phi_coord")) d.phi_coord = dict["phi_coord"].cast<double>();
            return d;
        }, "Create from Python dictionary")
        
        // Numpy array conversion
        .def("to_numpy", [](const MetricData &d) {
            // Return as structured numpy array
            py::array_t<double> arr = py::array_t<double>({1}, {
                sizeof(MetricData)
            });
            auto buf = arr.request();
            std::memcpy(buf.ptr, &d, sizeof(MetricData));
            return arr;
        }, "Convert to numpy array (advanced)");
}

// ============================================================================
// Core module bindings
// ============================================================================

void bind_core_module(py::module &m) {
    py::module core = m.def_submodule("core", "Core quasilocal functions");
    
    // Axisymmetric angular momentum
    core.def("compute_axisymmetric_j", 
        py::overload_cast<const std::vector<double>&,
                         const std::vector<double>&,
                         const std::vector<double>&>(
            &qlt::compute_axisymmetric_j),
        py::arg("alpha"), py::arg("f"), py::arg("r"),
        "Compute axisymmetric angular momentum J(r) = (1/3) α r² f")
    
    .def("compute_axisymmetric_j", 
        [](const py::array_t<double> &alpha,
           const py::array_t<double> &f,
           const py::array_t<double> &r) {
            // Numpy array version
            auto alpha_buf = alpha.request();
            auto f_buf = f.request();
            auto r_buf = r.request();
            
            if (alpha_buf.size != f_buf.size || f_buf.size != r_buf.size) {
                throw py::value_error("All arrays must have same size");
            }
            
            std::vector<double> alpha_vec(alpha_buf.size);
            std::vector<double> f_vec(f_buf.size);
            std::vector<double> r_vec(r_buf.size);
            
            std::memcpy(alpha_vec.data(), alpha_buf.ptr, alpha_buf.size * sizeof(double));
            std::memcpy(f_vec.data(), f_buf.ptr, f_buf.size * sizeof(double));
            std::memcpy(r_vec.data(), r_buf.ptr, r_buf.size * sizeof(double));
            
            auto result = qlt::compute_axisymmetric_j(alpha_vec, f_vec, r_vec);
            
            // Convert back to numpy array
            py::array_t<double> result_arr(result.size());
            auto result_buf = result_arr.request();
            std::memcpy(result_buf.ptr, result.data(), result.size() * sizeof(double));
            return result_arr;
        },
        py::arg("alpha").noconvert(),
        py::arg("f").noconvert(),
        py::arg("r").noconvert(),
        "Numpy array version of compute_axisymmetric_j")
    
    // Clebsch-Komar angular momentum
    .def("compute_clebsch_komar_j", &qlt::compute_clebsch_komar_j,
        py::arg("surface_data"), py::arg("volume_data"), py::arg("kappa") = 2.0,
        "Compute Clebsch-Komar angular momentum")
    
    // Quasilocal mass
    .def("compute_quasilocal_mass", &qlt::compute_quasilocal_mass,
        py::arg("surface_data"), py::arg("use_reference") = true,
        "Compute quasilocal mass")
    
    // Surgical flux
    .def("compute_surgical_flux", 
        py::overload_cast<const std::vector<MetricData>&,
                         const std::vector<MetricData>&,
                         bool>(&qlt::compute_surgical_flux),
        py::arg("inner_data"), py::arg("outer_data"), py::arg("axisymmetric") = false,
        "Compute surgical flux ΔJ across excision boundary")
    
    // Vorticity bounds
    .def("compute_vorticity_bound", &qlt::compute_vorticity_bound,
        py::arg("data"),
        "Compute vorticity bound: ||V ∧ dV|| ≤ ||α|| ||dβ||²")
    
    .def("check_vorticity_bound", &qlt::check_vorticity_bound,
        py::arg("data"), py::arg("tolerance") = 1e-10,
        "Check if vorticity bound is satisfied")
    
    // Zero angular momentum condition
    .def("check_zero_angular_momentum", &qlt::check_zero_angular_momentum,
        py::arg("data"), py::arg("tolerance") = 1e-10,
        "Check zero angular momentum condition: α≡0 or dβ≡0 ⇒ J=0")
    
    // Helicity
    .def("compute_helicity", &qlt::compute_helicity,
        py::arg("volume_data"),
        "Compute helicity ∫ V ∧ dV")
    
    // Integrated angular momentum bound
    .def("compute_integrated_angular_momentum_bound",
        py::overload_cast<const std::vector<MetricData>&, double, double>(
            &qlt::compute_integrated_angular_momentum_bound),
        py::arg("volume_data"), py::arg("adm_mass"), 
        py::arg("mass_constant_factor") = 1.0,
        "Compute integrated angular momentum bound")
    
    // ADM mass approximation
    .def("compute_adm_mass_approximation", &qlt::compute_adm_mass_approximation,
        py::arg("surface_data"),
        "Approximate ADM mass from asymptotic data")
    
    // Mass aspect
    .def("compute_mass_aspect", &qlt::compute_mass_aspect,
        py::arg("surface_data_at_radii"),
        "Compute mass aspect function M(r)");
}

// ============================================================================
// Grid module bindings
// ============================================================================

void bind_grid_module(py::module &m) {
    py::module grid = m.def_submodule("grid", "Grid data I/O and utilities");
    
    // DataFormat enum
    py::enum_<grid::DataFormat>(grid, "DataFormat")
        .value("SPECTRE_HDF5", grid::DataFormat::SPECTRE_HDF5)
        .value("CACTUS_HDF5", grid::DataFormat::CACTUS_HDF5)
        .value("CACTUS_ASCII", grid::DataFormat::CACTUS_ASCII)
        .value("SXS_CATALOG", grid::DataFormat::SXS_CATALOG)
        .value("CUSTOM_CSV", grid::DataFormat::CUSTOM_CSV)
        .value("UNKNOWN", grid::DataFormat::UNKNOWN)
        .export_values();
    
    // GridData class
    py::class_<grid::GridData>(grid, "GridData", "3D grid data structure")
        .def(py::init<>())
        .def_readwrite("data", &grid::GridData::data)
        .def_readwrite("x_coords", &grid::GridData::x_coords)
        .def_readwrite("y_coords", &grid::GridData::y_coords)
        .def_readwrite("z_coords", &grid::GridData::z_coords)
        .def_readwrite("coordinate_system", &grid::GridData::coordinate_system)
        .def_readwrite("time", &grid::GridData::time)
        .def_readwrite("iteration", &grid::GridData::iteration)
        .def_readwrite("format", &grid::GridData::format)
        .def_readwrite("filename", &grid::GridData::filename)
        .def("valid", &grid::GridData::valid, "Check if grid is valid")
        .def("dimensions", &grid::GridData::dimensions, "Get grid dimensions")
        .def("__repr__", [](const grid::GridData &g) {
            std::ostringstream oss;
            auto dims = g.dimensions();
            oss << "GridData(dimensions=[" << dims[0] << ", " 
                << dims[1] << ", " << dims[2] << "], "
                << "time=" << g.time << ")";
            return oss.str();
        });
    
    // MetricReader class
    py::class_<grid::MetricReader>(grid, "MetricReader", 
        "Read numerical relativity data from files")
        .def(py::init<>())
        
        .def_static("detect_format", &grid::MetricReader::detect_format,
            py::arg("filename"), "Detect data format from filename")
        
        .def("read_file", &grid::MetricReader::read_file,
            py::arg("filename"), py::arg("format") = grid::DataFormat::UNKNOWN,
            "Read data from file")
        
        .def("read_timeslice", &grid::MetricReader::read_timeslice,
            py::arg("filename"), py::arg("time"), 
            py::arg("format") = grid::DataFormat::SPECTRE_HDF5,
            "Read specific time slice")
        
        .def("extract_spherical_shell", &grid::MetricReader::extract_spherical_shell,
            py::arg("cartesian_data"), py::arg("radius"),
            py::arg("num_theta") = 64, py::arg("num_phi") = 128,
            "Extract spherical shell from Cartesian data")
        
        .def("get_metadata", &grid::MetricReader::get_metadata,
            py::arg("filename"), "Get metadata from file")
        
        .def("check_required_fields", &grid::MetricReader::check_required_fields,
            py::arg("filename"), py::arg("fields"),
            "Check if file contains required fields");
    
    // CoordinateSystems class
    py::class_<grid::CoordinateSystems>(grid, "CoordinateSystems",
        "Coordinate transformation utilities")
        .def_static("cartesian_to_spherical", 
            py::overload_cast<double, double, double>(
                &grid::CoordinateSystems::cartesian_to_spherical),
            py::arg("x"), py::arg("y"), py::arg("z"),
            "Convert Cartesian to spherical coordinates")
        
        .def_static("spherical_to_cartesian",
            py::overload_cast<double, double, double>(
                &grid::CoordinateSystems::spherical_to_cartesian),
            py::arg("r"), py::arg("theta"), py::arg("phi"),
            "Convert spherical to Cartesian coordinates")
        
        .def_static("generate_spherical_grid",
            &grid::CoordinateSystems::generate_spherical_grid,
            py::arg("r_min"), py::arg("r_max"), py::arg("num_r"),
            py::arg("num_theta"), py::arg("num_phi"),
            "Generate spherical grid coordinates")
        
        .def_static("validate_coordinates",
            &grid::CoordinateSystems::validate_coordinates,
            py::arg("coords"), py::arg("system"),
            "Validate coordinates for given system");
    
    // KillingVectorType enum
    py::enum_<grid::KillingVectorType>(grid, "KillingVectorType")
        .value("ROTATIONAL_Z", grid::KillingVectorType::ROTATIONAL_Z)
        .value("ROTATIONAL_X", grid::KillingVectorType::ROTATIONAL_X)
        .value("ROTATIONAL_Y", grid::KillingVectorType::ROTATIONAL_Y)
        .value("TRANSLATIONAL_X", grid::KillingVectorType::TRANSLATIONAL_X)
        .value("TRANSLATIONAL_Y", grid::KillingVectorType::TRANSLATIONAL_Y)
        .value("TRANSLATIONAL_Z", grid::KillingVectorType::TRANSLATIONAL_Z)
        .value("BOOST_X", grid::KillingVectorType::BOOST_X)
        .value("BOOST_Y", grid::KillingVectorType::BOOST_Y)
        .value("BOOST_Z", grid::KillingVectorType::BOOST_Z)
        .value("CUSTOM", grid::KillingVectorType::CUSTOM)
        .value("APPROXIMATE", grid::KillingVectorType::APPROXIMATE)
        .export_values();
    
    // KillingVectorGenerator class
    py::class_<grid::KillingVectorGenerator>(grid, "KillingVectorGenerator",
        "Generate Killing vector fields")
        .def(py::init<>())
        
        .def("generate_killing_field", &grid::KillingVectorGenerator::generate_killing_field,
            py::arg("grid_data"), py::arg("type"), 
            py::arg("parameters") = std::array<double, 3>{0.0, 0.0, 0.0},
            "Generate Killing vector field of specified type")
        
        .def("generate_rotational_killing_vector",
            &grid::KillingVectorGenerator::generate_rotational_killing_vector,
            py::arg("grid_data"),
            py::arg("axis") = std::array<double, 3>{0.0, 0.0, 1.0},
            py::arg("angular_velocity") = 1.0,
            "Generate rotational Killing vector about specified axis")
        
        .def("check_killing_equation", 
            [](grid::KillingVectorGenerator &self,
               const std::vector<MetricData> &data_with_xi) {
                double max_violation, rms_violation;
                auto result = self.check_killing_equation(
                    data_with_xi, max_violation, rms_violation);
                return py::make_tuple(result.first, result.second,
                                     max_violation, rms_violation);
            },
            py::arg("data_with_xi"),
            "Check Killing equation violation")
        
        .def("save_to_file", &grid::KillingVectorGenerator::save_to_file,
            py::arg("data_with_xi"), py::arg("filename"),
            "Save Killing vector to file");
    
    // SurfaceIntegral class
    py::class_<grid::SurfaceIntegral>(grid, "SurfaceIntegral",
        "Surface integration utilities")
        .def_static("compute_surface_area", &grid::SurfaceIntegral::compute_surface_area,
            py::arg("surface_data"), "Compute area of 2-surface")
        
        .def_static("compute_komar_surface_integral",
            &grid::SurfaceIntegral::compute_komar_surface_integral,
            py::arg("surface_data"), "Compute Komar surface integral")
        
        .def_static("compute_surgical_flux_integral",
            &grid::SurfaceIntegral::compute_surgical_flux_integral,
            py::arg("surface_data"),
            "Compute surgical flux surface integral");
    
    // VolumeIntegral class
    py::class_<grid::VolumeIntegral>(grid, "VolumeIntegral",
        "Volume integration utilities")
        .def_static("compute_volume", &grid::VolumeIntegral::compute_volume,
            py::arg("volume_data"), "Compute volume of region")
        
        .def_static("compute_clebsch_komar_volume_integral",
            &grid::VolumeIntegral::compute_clebsch_komar_volume_integral,
            py::arg("volume_data"),
            "Compute Clebsch-Komar volume integral")
        
        .def_static("compute_vorticity_bound_integral",
            &grid::VolumeIntegral::compute_vorticity_bound_integral,
            py::arg("volume_data"),
            "Compute vorticity bound volume integral");
    
    // Helper function to create MetricData from numpy arrays
    grid.def("create_metric_data_from_arrays", 
        [](const py::array_t<double> &h_ij,
           const py::array_t<double> &K_ij,
           const py::array_t<double> &coords,
           const py::array_t<double> &xi,
           const py::array_t<double> &dxi_dx,
           double alpha, double phi, double beta,
           const py::array_t<double> &dalpha_dx,
           const py::array_t<double> &dbeta_dx,
           const py::array_t<double> &dphi_dx) {
            
            // Convert numpy arrays to C++ pointers
            double* h_ij_ptr = nullptr;
            double* K_ij_ptr = nullptr;
            double* coords_ptr = nullptr;
            double* xi_ptr = nullptr;
            double* dxi_dx_ptr = nullptr;
            double* dalpha_dx_ptr = nullptr;
            double* dbeta_dx_ptr = nullptr;
            double* dphi_dx_ptr = nullptr;
            
            if (!h_ij.is_none()) {
                auto h_ij_buf = h_ij.request();
                if (h_ij_buf.size == 9) h_ij_ptr = static_cast<double*>(h_ij_buf.ptr);
            }
            
            if (!K_ij.is_none()) {
                auto K_ij_buf = K_ij.request();
                if (K_ij_buf.size == 9) K_ij_ptr = static_cast<double*>(K_ij_buf.ptr);
            }
            
            if (!coords.is_none()) {
                auto coords_buf = coords.request();
                if (coords_buf.size == 3) coords_ptr = static_cast<double*>(coords_buf.ptr);
            }
            
            // Similar for other arrays...
            
            return grid::create_metric_data_from_arrays(
                h_ij_ptr, K_ij_ptr, coords_ptr,
                xi_ptr, dxi_dx_ptr,
                alpha, phi, beta,
                dalpha_dx_ptr, dbeta_dx_ptr, dphi_dx_ptr);
        },
        py::arg("h_ij") = py::none(), py::arg("K_ij") = py::none(),
        py::arg("coords") = py::none(), py::arg("xi") = py::none(),
        py::arg("dxi_dx") = py::none(), py::arg("alpha") = 0.0,
        py::arg("phi") = 0.0, py::arg("beta") = 0.0,
        py::arg("dalpha_dx") = py::none(), py::arg("dbeta_dx") = py::none(),
        py::arg("dphi_dx") = py::none(),
        "Create MetricData from numpy arrays");
}

} // namespace python
} // namespace qlt

// ============================================================================
// Main module definition
// ============================================================================

PYBIND11_MODULE(qlt, m) {
    m.doc() = R"pbdoc(
        QuasiLocal Toolkit (QLT)
        -------------------------
        
        Mathematically accurate implementation of quasilocal conserved quantities
        from "Geometric Vorticity and Quasilocal Angular Momentum Conservation"
        by Rayan D. Peru.
        
        Provides tools for computing angular momentum, mass, and other conserved
        quantities in numerical relativity without horizon finding.
    )pbdoc";
    
    // Add version info
    m.attr("__version__") = "0.1.0";
    m.attr("__author__") = "Rayan D. Peru";
    
    // Add constants
    m.attr("M_PI") = M_PI;
    
    // Bind all modules
    qlt::python::bind_metric_data(m);
    qlt::python::bind_core_module(m);
    qlt::python::bind_grid_module(m);
    
    // Add convenience functions
    m.def("test_installation", []() {
        py::dict result;
        result["status"] = "OK";
        result["version"] = "0.1.0";
        result["message"] = "QLT Python bindings installed successfully";
        return result;
    }, "Test installation");
    
    m.def("about", []() {
        std::ostringstream oss;
        oss << "QuasiLocal Toolkit (QLT) v0.1.0\n"
            << "Mathematical implementation of quasilocal conserved quantities\n"
            << "Based on: Geometric Vorticity and Quasilocal Angular Momentum Conservation\n"
            << "Author: Rayan D. Peru\n"
            << "License: MIT";
        return oss.str();
    }, "About QLT");
}