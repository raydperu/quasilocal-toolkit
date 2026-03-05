#ifndef QLT_METRIC_READER_H
#define QLT_METRIC_READER_H

#include "../core/quasilocal_quantities.h"
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace qlt {
namespace grid {

// Forward declarations for PIMPL pattern
struct HDF5Impl;
struct SpECDataImpl;
struct CactusDataImpl;

/**
 * @brief Supported numerical relativity data formats
 */
enum class DataFormat {
    SPECTRE_HDF5,   // SpEC/SXS HDF5 format
    CACTUS_HDF5,    // Einstein Toolkit/Cactus HDF5
    CACTUS_ASCII,   // Cactus ASCII output
    SXS_CATALOG,    // SXS waveform catalog
    CUSTOM_CSV,     // Custom CSV format
    UNKNOWN
};

/**
 * @brief Grid data structure for 3D numerical relativity data
 */
struct GridData {
    std::vector<std::vector<std::vector<MetricData>>> data;  // 3D grid [i][j][k]
    
    // Grid dimensions and coordinates
    std::vector<double> x_coords;  // x_i values
    std::vector<double> y_coords;  // y_j values
    std::vector<double> z_coords;  // z_k values
    
    // Coordinate system information
    std::string coordinate_system;  // "Cartesian", "Spherical", etc.
    
    // Time slice information
    double time;                    // Simulation time
    int iteration;                  // Iteration number
    
    // Metadata
    DataFormat format;
    std::string filename;
    
    // Check if grid is valid
    bool valid() const {
        return !data.empty() && !data[0].empty() && !data[0][0].empty();
    }
    
    // Get grid dimensions
    std::array<size_t, 3> dimensions() const {
        if (!valid()) return {0, 0, 0};
        return {
            data.size(),
            data[0].size(),
            data[0][0].size()
        };
    }
};

/**
 * @brief Main class for reading numerical relativity data
 */
class MetricReader {
public:
    MetricReader();
    ~MetricReader();
    
    // Disable copy (use PIMPL)
    MetricReader(const MetricReader&) = delete;
    MetricReader& operator=(const MetricReader&) = delete;
    
    // Move operations
    MetricReader(MetricReader&&) noexcept;
    MetricReader& operator=(MetricReader&&) noexcept;
    
    /**
     * @brief Detect data format from filename
     */
    static DataFormat detect_format(const std::string& filename);
    
    /**
     * @brief Read data from file
     * @param filename Input file
     * @param format Optional format specification
     * @return GridData containing the read data
     */
    GridData read_file(const std::string& filename, 
                      DataFormat format = DataFormat::UNKNOWN);
    
    /**
     * @brief Read specific time slice from HDF5 file
     */
    GridData read_timeslice(const std::string& filename,
                           double time,
                           DataFormat format = DataFormat::SPECTRE_HDF5);
    
    /**
     * @brief Read multiple time slices
     */
    std::vector<GridData> read_timeseries(const std::string& filename,
                                         const std::vector<double>& times);
    
    /**
     * @brief Extract spherical shell from Cartesian data
     * @param cartesian_data Input Cartesian grid
     * @param radius Desired radius
     * @param num_theta Number of θ points
     * @param num_phi Number of φ points
     * @return Data interpolated onto sphere
     */
    std::vector<MetricData> extract_spherical_shell(
        const GridData& cartesian_data,
        double radius,
        size_t num_theta = 64,
        size_t num_phi = 128);
    
    /**
     * @brief Extract cylindrical surface
     */
    std::vector<MetricData> extract_cylindrical_surface(
        const GridData& cartesian_data,
        double radius,
        double z_min,
        double z_max,
        size_t num_phi = 128,
        size_t num_z = 64);
    
    /**
     * @brief Get metadata from file
     */
    std::map<std::string, std::string> get_metadata(
        const std::string& filename);
    
    /**
     * @brief Check if file contains required fields
     */
    bool check_required_fields(const std::string& filename,
                              const std::vector<std::string>& fields);
    
private:
    std::unique_ptr<HDF5Impl> hdf5_impl_;
    std::unique_ptr<SpECDataImpl> spec_impl_;
    std::unique_ptr<CactusDataImpl> cactus_impl_;
    
    // Format-specific readers
    GridData read_spectre_hdf5(const std::string& filename);
    GridData read_cactus_hdf5(const std::string& filename);
    GridData read_cactus_ascii(const std::string& filename);
    GridData read_sxs_catalog(const std::string& filename);
    GridData read_custom_csv(const std::string& filename);
};

/**
 * @brief Helper function to create MetricData from raw arrays
 */
MetricData create_metric_data_from_arrays(
    const double* h_ij,      // 9 elements (3x3)
    const double* K_ij,      // 9 elements (3x3)
    const double* coords,    // 3 elements (x,y,z)
    const double* xi = nullptr,      // 3 elements (optional)
    const double* dxi_dx = nullptr,  // 9 elements (optional)
    double alpha = 0.0,
    double phi = 0.0,
    double beta = 0.0,
    const double* dalpha_dx = nullptr,
    const double* dbeta_dx = nullptr,
    const double* dphi_dx = nullptr);

} // namespace grid
} // namespace qlt

#endif // QLT_METRIC_READER_H