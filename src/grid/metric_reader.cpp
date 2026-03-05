#include "grid/metric_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <regex>

// HDF5 includes (conditional compilation)
#ifdef HAVE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

namespace qlt {
namespace grid {

// ============================================================================
// PIMPL Implementation Structures
// ============================================================================

struct HDF5Impl {
#ifdef HAVE_HDF5
    hid_t file_id = -1;
    hid_t group_id = -1;
    hid_t dataset_id = -1;
#endif
    
    ~HDF5Impl() {
#ifdef HAVE_HDF5
        if (dataset_id >= 0) H5Dclose(dataset_id);
        if (group_id >= 0) H5Gclose(group_id);
        if (file_id >= 0) H5Fclose(file_id);
#endif
    }
};

struct SpECDataImpl {
    // SpEC-specific data structures
    std::map<std::string, std::vector<double>> datasets;
    std::map<std::string, std::vector<int>> attributes;
};

struct CactusDataImpl {
    // Cactus-specific data structures
    std::vector<std::string> variable_names;
    std::map<std::string, std::vector<double>> data_arrays;
};

// ============================================================================
// MetricReader Implementation
// ============================================================================

MetricReader::MetricReader() 
    : hdf5_impl_(std::make_unique<HDF5Impl>())
    , spec_impl_(std::make_unique<SpECDataImpl>())
    , cactus_impl_(std::make_unique<CactusDataImpl>()) {
}

MetricReader::~MetricReader() = default;

MetricReader::MetricReader(MetricReader&&) noexcept = default;
MetricReader& MetricReader::operator=(MetricReader&&) noexcept = default;

DataFormat MetricReader::detect_format(const std::string& filename) {
    std::string extension;
    size_t dot_pos = filename.find_last_of('.');
    if (dot_pos != std::string::npos) {
        extension = filename.substr(dot_pos + 1);
        std::transform(extension.begin(), extension.end(), 
                      extension.begin(), ::tolower);
    }
    
    // Check file signature for HDF5
    std::ifstream file(filename, std::ios::binary);
    if (file) {
        char signature[4];
        file.read(signature, 4);
        if (std::memcmp(signature, "\211HDF\r\n\032\n", 8) == 0) {
            // Check if it's SpEC format by looking for typical groups
            return DataFormat::SPECTRE_HDF5;
        }
    }
    
    // Check extension
    if (extension == "h5" || extension == "hdf5") {
        return DataFormat::SPECTRE_HDF5;
    } else if (extension == "asc" || extension == "dat") {
        return DataFormat::CACTUS_ASCII;
    } else if (extension == "csv") {
        return DataFormat::CUSTOM_CSV;
    }
    
    return DataFormat::UNKNOWN;
}

GridData MetricReader::read_file(const std::string& filename, 
                                DataFormat format) {
    
    if (format == DataFormat::UNKNOWN) {
        format = detect_format(filename);
    }
    
    std::cout << "[INFO] Reading file: " << filename 
              << " (format: " << static_cast<int>(format) << ")" << std::endl;
    
    switch (format) {
        case DataFormat::SPECTRE_HDF5:
            return read_spectre_hdf5(filename);
        case DataFormat::CACTUS_HDF5:
            return read_cactus_hdf5(filename);
        case DataFormat::CACTUS_ASCII:
            return read_cactus_ascii(filename);
        case DataFormat::SXS_CATALOG:
            return read_sxs_catalog(filename);
        case DataFormat::CUSTOM_CSV:
            return read_custom_csv(filename);
        default:
            throw std::runtime_error("Unsupported data format for file: " + filename);
    }
}

GridData MetricReader::read_spectre_hdf5(const std::string& filename) {
    GridData grid_data;
    grid_data.format = DataFormat::SPECTRE_HDF5;
    grid_data.filename = filename;
    
#ifdef HAVE_HDF5
    std::cout << "[INFO] Reading SpEC HDF5 file: " << filename << std::endl;
    
    // Open HDF5 file
    hdf5_impl_->file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (hdf5_impl_->file_id < 0) {
        throw std::runtime_error("Failed to open HDF5 file: " + filename);
    }
    
    // Read metadata
    herr_t status;
    H5Eset_auto(NULL, NULL); // Turn off error printing
    
    // Try to read different group structures (SpEC evolution output)
    const char* group_paths[] = {
        "/",
        "/AhA.dir",
        "/AhB.dir",
        "/VolumeData",
        "/ApparentHorizons"
    };
    
    for (const char* path : group_paths) {
        if (H5Lexists(hdf5_impl_->file_id, path, H5P_DEFAULT) > 0) {
            hdf5_impl_->group_id = H5Gopen(hdf5_impl_->file_id, path, H5P_DEFAULT);
            if (hdf5_impl_->group_id >= 0) {
                std::cout << "[INFO] Found group: " << path << std::endl;
                break;
            }
        }
    }
    
    if (hdf5_impl_->group_id < 0) {
        // Try to find any group
        hsize_t num_objs;
        H5Gget_num_objs(hdf5_impl_->file_id, &num_objs);
        for (hsize_t i = 0; i < num_objs; i++) {
            char name[256];
            H5Gget_objname_by_idx(hdf5_impl_->file_id, i, name, 256);
            if (H5Gget_objtype_by_idx(hdf5_impl_->file_id, i) == H5G_GROUP) {
                hdf5_impl_->group_id = H5Gopen(hdf5_impl_->file_id, name, H5P_DEFAULT);
                if (hdf5_impl_->group_id >= 0) {
                    std::cout << "[INFO] Using group: " << name << std::endl;
                    break;
                }
            }
        }
    }
    
    if (hdf5_impl_->group_id < 0) {
        H5Fclose(hdf5_impl_->file_id);
        throw std::runtime_error("No valid group found in HDF5 file");
    }
    
    // Read datasets (simplified - in practice would read specific variables)
    // gxx, gxy, gxz, gyy, gyz, gzz, Kxx, Kxy, etc.
    
    // For now, create dummy data for testing
    std::cout << "[WARNING] Using dummy data for SpEC HDF5 reader" << std::endl;
    
    // Create a simple 8x8x8 grid for testing
    size_t nx = 8, ny = 8, nz = 8;
    grid_data.x_coords.resize(nx);
    grid_data.y_coords.resize(ny);
    grid_data.z_coords.resize(nz);
    
    // Fill coordinates
    for (size_t i = 0; i < nx; i++) grid_data.x_coords[i] = -1.0 + 2.0 * i / (nx - 1);
    for (size_t j = 0; j < ny; j++) grid_data.y_coords[j] = -1.0 + 2.0 * j / (ny - 1);
    for (size_t k = 0; k < nz; k++) grid_data.z_coords[k] = -1.0 + 2.0 * k / (nz - 1);
    
    // Create grid data structure
    grid_data.data.resize(nx);
    for (size_t i = 0; i < nx; i++) {
        grid_data.data[i].resize(ny);
        for (size_t j = 0; j < ny; j++) {
            grid_data.data[i][j].resize(nz);
            for (size_t k = 0; k < nz; k++) {
                MetricData& md = grid_data.data[i][j][k];
                
                // Simple flat metric for testing
                md.h_ij = Matrix3d::Identity();
                md.K_ij = Matrix3d::Zero();
                
                // Position
                md.r = std::sqrt(
                    grid_data.x_coords[i] * grid_data.x_coords[i] +
                    grid_data.y_coords[j] * grid_data.y_coords[j] +
                    grid_data.z_coords[k] * grid_data.z_coords[k]
                );
                
                // Spherical coordinates
                if (md.r > 1e-10) {
                    md.theta = std::acos(grid_data.z_coords[k] / md.r);
                    md.phi_coord = std::atan2(grid_data.y_coords[j], grid_data.x_coords[i]);
                } else {
                    md.theta = 0.0;
                    md.phi_coord = 0.0;
                }
                
                // Default Clebsch potentials (for Kerr test)
                md.alpha = 0.0;
                md.phi = 0.0;
                md.beta = md.phi_coord; // β = φ in flat space
                
                // Derivatives (zero for flat space)
                md.dalpha_dx = Vector3d::Zero();
                md.dbeta_dx = Vector3d::Zero();
                md.dphi_dx = Vector3d::Zero();
                
                // Killing vector (rotational about z-axis)
                md.xi = Vector3d(-grid_data.y_coords[j], grid_data.x_coords[i], 0.0);
                md.dxi_dx = Matrix3d::Zero();
                md.dxi_dx(0,1) = -1.0;  // ∂ξ^x/∂y = -1
                md.dxi_dx(1,0) = 1.0;   // ∂ξ^y/∂x = 1
                
                // Metric determinant
                md.sqrt_det_h = 1.0;
            }
        }
    }
    
    grid_data.coordinate_system = "Cartesian";
    grid_data.time = 0.0;
    grid_data.iteration = 0;
    
    // Cleanup
    H5Gclose(hdf5_impl_->group_id);
    H5Fclose(hdf5_impl_->file_id);
    
#else
    throw std::runtime_error("HDF5 support not compiled in");
#endif
    
    return grid_data;
}

GridData MetricReader::read_cactus_hdf5(const std::string& filename) {
    GridData grid_data;
    grid_data.format = DataFormat::CACTUS_HDF5;
    grid_data.filename = filename;
    
    std::cout << "[WARNING] Cactus HDF5 reader not fully implemented" << std::endl;
    std::cout << "[INFO] Creating dummy data for testing" << std::endl;
    
    // Create dummy data similar to SpEC reader
    return read_spectre_hdf5(filename); // Reuse for now
}

GridData MetricReader::read_cactus_ascii(const std::string& filename) {
    GridData grid_data;
    grid_data.format = DataFormat::CACTUS_ASCII;
    grid_data.filename = filename;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::cout << "[INFO] Reading Cactus ASCII file: " << filename << std::endl;
    
    // Parse header
    std::string line;
    std::vector<std::string> column_names;
    size_t num_points = 0;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            // Comment line, check for metadata
            if (line.find("num_points") != std::string::npos) {
                std::istringstream iss(line);
                std::string token;
                while (iss >> token) {
                    if (token.find("=") != std::string::npos) {
                        size_t eq_pos = token.find('=');
                        std::string key = token.substr(0, eq_pos);
                        std::string value = token.substr(eq_pos + 1);
                        if (key == "num_points") {
                            num_points = std::stoul(value);
                        }
                    }
                }
            }
            continue;
        }
        
        // First non-comment line is header
        std::istringstream iss(line);
        std::string token;
        while (iss >> token) {
            column_names.push_back(token);
        }
        break;
    }
    
    std::cout << "[INFO] Found " << column_names.size() << " columns" << std::endl;
    std::cout << "[INFO] num_points = " << num_points << std::endl;
    
    // For now, create dummy data
    // In full implementation, would parse the actual data
    
    return grid_data;
}

GridData MetricReader::read_sxs_catalog(const std::string& filename) {
    GridData grid_data;
    grid_data.format = DataFormat::SXS_CATALOG;
    grid_data.filename = filename;
    
    std::cout << "[INFO] SXS catalog reader not implemented" << std::endl;
    // Would connect to SXS catalog API
    // https://data.black-holes.org
    
    return grid_data;
}

GridData MetricReader::read_custom_csv(const std::string& filename) {
    GridData grid_data;
    grid_data.format = DataFormat::CUSTOM_CSV;
    grid_data.filename = filename;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open CSV file: " + filename);
    }
    
    std::cout << "[INFO] Reading custom CSV file: " << filename << std::endl;
    
    // Parse CSV
    std::string line;
    std::vector<std::string> headers;
    std::vector<std::vector<double>> columns;
    
    // Read header
    if (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        while (std::getline(iss, token, ',')) {
            headers.push_back(token);
            columns.emplace_back();
        }
    }
    
    // Read data
    size_t row = 0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        size_t col = 0;
        
        while (std::getline(iss, token, ',') && col < headers.size()) {
            try {
                double value = std::stod(token);
                columns[col].push_back(value);
            } catch (const std::exception& e) {
                std::cerr << "[WARNING] Could not parse value at row " << row 
                          << ", col " << col << ": " << token << std::endl;
                columns[col].push_back(0.0);
            }
            col++;
        }
        row++;
    }
    
    std::cout << "[INFO] Read " << row << " rows with " << headers.size() << " columns" << std::endl;
    
    // Convert to GridData (simplified 1D array for now)
    // In full implementation, would need to know grid structure
    
    return grid_data;
}

std::vector<MetricData> MetricReader::extract_spherical_shell(
    const GridData& cartesian_data,
    double radius,
    size_t num_theta,
    size_t num_phi) {
    
    if (!cartesian_data.valid()) {
        throw std::invalid_argument("Invalid Cartesian data");
    }
    
    if (cartesian_data.coordinate_system != "Cartesian") {
        std::cerr << "[WARNING] Input data not in Cartesian coordinates" << std::endl;
    }
    
    std::vector<MetricData> shell_data;
    shell_data.reserve(num_theta * num_phi);
    
    // Create spherical grid
    for (size_t i = 0; i < num_theta; i++) {
        double theta = M_PI * i / (num_theta - 1);
        if (i == 0) theta = 1e-10;  // Avoid pole
        if (i == num_theta - 1) theta = M_PI - 1e-10;
        
        for (size_t j = 0; j < num_phi; j++) {
            double phi = 2.0 * M_PI * j / num_phi;
            
            // Convert to Cartesian for interpolation
            double x = radius * sin(theta) * cos(phi);
            double y = radius * sin(theta) * sin(phi);
            double z = radius * cos(theta);
            
            // Find nearest grid point (simplified - would use trilinear interpolation)
            size_t ix = std::min(cartesian_data.x_coords.size() - 1,
                               static_cast<size_t>((x - cartesian_data.x_coords[0]) / 
                               (cartesian_data.x_coords[1] - cartesian_data.x_coords[0])));
            size_t iy = std::min(cartesian_data.y_coords.size() - 1,
                               static_cast<size_t>((y - cartesian_data.y_coords[0]) / 
                               (cartesian_data.y_coords[1] - cartesian_data.y_coords[0])));
            size_t iz = std::min(cartesian_data.z_coords.size() - 1,
                               static_cast<size_t>((z - cartesian_data.z_coords[0]) / 
                               (cartesian_data.z_coords[1] - cartesian_data.z_coords[0])));
            
            if (ix < cartesian_data.data.size() &&
                iy < cartesian_data.data[ix].size() &&
                iz < cartesian_data.data[ix][iy].size()) {
                
                MetricData data = cartesian_data.data[ix][iy][iz];
                
                // Update position for spherical shell
                data.r = radius;
                data.theta = theta;
                data.phi_coord = phi;
                
                // Compute normal vector (radial outward)
                data.normal = Vector3d(sin(theta) * cos(phi),
                                      sin(theta) * sin(phi),
                                      cos(theta));
                
                shell_data.push_back(data);
            }
        }
    }
    
    std::cout << "[INFO] Extracted spherical shell with radius " << radius
              << ": " << shell_data.size() << " points" << std::endl;
    
    return shell_data;
}

std::vector<MetricData> MetricReader::extract_cylindrical_surface(
    const GridData& cartesian_data,
    double radius,
    double z_min,
    double z_max,
    size_t num_phi,
    size_t num_z) {
    
    std::vector<MetricData> surface_data;
    surface_data.reserve(num_phi * num_z);
    
    // Create cylindrical grid
    for (size_t i = 0; i < num_phi; i++) {
        double phi = 2.0 * M_PI * i / num_phi;
        
        for (size_t j = 0; j < num_z; j++) {
            double z = z_min + (z_max - z_min) * j / (num_z - 1);
            
            // Convert to Cartesian
            double x = radius * cos(phi);
            double y = radius * sin(phi);
            
            // Find nearest grid point (simplified)
            // ... similar to spherical shell extraction
            
            // For now, create dummy data
            MetricData data;
            data.r = radius;
            data.phi_coord = phi;
            
            // Cylindrical normal (radial outward)
            data.normal = Vector3d(cos(phi), sin(phi), 0.0);
            
            surface_data.push_back(data);
        }
    }
    
    return surface_data;
}

std::map<std::string, std::string> MetricReader::get_metadata(
    const std::string& filename) {
    
    std::map<std::string, std::string> metadata;
    DataFormat format = detect_format(filename);
    
    if (format == DataFormat::SPECTRE_HDF5) {
#ifdef HAVE_HDF5
        // Read HDF5 attributes
#endif
        metadata["format"] = "SpEC_HDF5";
        metadata["version"] = "1.0";
    } else if (format == DataFormat::CACTUS_ASCII) {
        metadata["format"] = "Cactus_ASCII";
    }
    
    // Add file info
    metadata["filename"] = filename;
    
    return metadata;
}

bool MetricReader::check_required_fields(const std::string& filename,
                                        const std::vector<std::string>& fields) {
    // Simplified check
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    
    // For HDF5, would check for datasets
    // For ASCII, check header
    
    return true;
}

MetricData create_metric_data_from_arrays(
    const double* h_ij,
    const double* K_ij,
    const double* coords,
    const double* xi,
    const double* dxi_dx,
    double alpha,
    double phi,
    double beta,
    const double* dalpha_dx,
    const double* dbeta_dx,
    const double* dphi_dx) {
    
    MetricData data;
    
    // Copy metric
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            data.h_ij(i, j) = h_ij[i * 3 + j];
            data.K_ij(i, j) = K_ij[i * 3 + j];
        }
    }
    
    // Copy coordinates
    if (coords) {
        // Compute spherical coordinates if needed
        double x = coords[0], y = coords[1], z = coords[2];
        data.r = std::sqrt(x * x + y * y + z * z);
        if (data.r > 1e-10) {
            data.theta = std::acos(z / data.r);
            data.phi_coord = std::atan2(y, x);
        } else {
            data.theta = 0.0;
            data.phi_coord = 0.0;
        }
    }
    
    // Copy Killing vector
    if (xi) {
        data.xi = Vector3d(xi[0], xi[1], xi[2]);
    }
    
    // Copy Killing vector derivative
    if (dxi_dx) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                data.dxi_dx(i, j) = dxi_dx[i * 3 + j];
            }
        }
    }
    
    // Clebsch potentials
    data.alpha = alpha;
    data.phi = phi;
    data.beta = beta;
    
    // Clebsch potential derivatives
    if (dalpha_dx) {
        data.dalpha_dx = Vector3d(dalpha_dx[0], dalpha_dx[1], dalpha_dx[2]);
    }
    if (dbeta_dx) {
        data.dbeta_dx = Vector3d(dbeta_dx[0], dbeta_dx[1], dbeta_dx[2]);
    }
    if (dphi_dx) {
        data.dphi_dx = Vector3d(dphi_dx[0], dphi_dx[1], dphi_dx[2]);
    }
    
    // Compute metric determinant
    data.sqrt_det_h = std::sqrt(data.h_ij.determinant());
    
    return data;
}

} // namespace grid
} // namespace qlt