# IsoHysplit

IsoHysplit is a Python package for calculating and analyzing HYSPLIT back trajectories with a focus on water isotope data integration. It provides tools for trajectory generation, analysis, and visualization including moisture source diagnostics and isotopic composition along air parcel trajectories.

## Acknowledgments

This code benefits significantly from the following projects:
- [splitr](https://github.com/rich-iannone/splitr)
- [HyTraj](https://github.com/pankajkarman/HyTraj) 
- [PySPLIT](https://github.com/mscross/pysplit)

## Requirements

- Python 3.x
- HYSPLIT installation
- Required Python packages:
  ```
  numpy
  pandas
  xarray
  matplotlib
  cartopy
  geopandas
  shapely
  joblib
  ```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/IsoHysplit.git
   cd IsoHysplit
   ```

2. Install required Python packages:
   ```bash
   pip install -r requirements.txt
   ```

3. Ensure HYSPLIT is installed and the executable path is correctly specified in your namelist file.

## Directory Structure

```
IsoHysplit/
├── exec/                  # HYSPLIT executables
├── working/              # Working directory for HYSPLIT
├── output/               # Output directory
├── Forcing/              # Meteorological data
├── isotope/              # Isotope data files
├── namelist/             # Namelist configuration files
└── src/                  # Source code
```

## Configuration

Configure your analysis by editing the namelist file (main.nml). Example configuration:

```fortran
&general
  basename=Singapore
  years=2005
  months=4
  hours = 0, 3, 6, 9, 12, 15, 18, 21
  altitudes = 1000
  location = 1.352, 103.820
  runtime = -240
  working_dir   = ../working
  storage_dir   = ../output/SG2
  meteo_dir     = ../Forcing/gdas1-new
  isotope_dir   = ../isotope
  hysplit       = ../exec
/
```

Key parameters include:
- `basename`: Base name for output files
- `years`, `months`, `hours`: Temporal parameters
- `altitudes`: Starting heights for trajectories (in meters)
- `location`: Starting coordinates (latitude, longitude)
- `runtime`: Length of back trajectories in hours (negative for back trajectories)
- Directory paths:
  - `working_dir`: HYSPLIT working directory
  - `storage_dir`: Output storage location
  - `meteo_dir`: Meteorological data location
  - `isotope_dir`: Isotope data location

## Usage

1. Prepare your namelist file:
   ```bash
   cp namelist/main.nml.template namelist/main.nml
   # Edit main.nml with your parameters
   ```

2. Run the main script:
   ```bash
   python IsoHysplit.py namelist/main.nml
   ```

## Available Analyses

### 1. Trajectory Generation
- Back trajectory calculation using HYSPLIT
- Multiple starting heights and times
- Configurable runtime and meteorological inputs

### 2. Basic Trajectory Analysis
- Humidity along trajectories
- Pressure levels
- Moisture flux calculations
- Moisture uptake diagnostics

### 3. Isotope Analysis
- Precipitable water δ18O and δD analysis
- Precipitation water δ18O and δD analysis
- Vapor water δ18O and δD analysis
- Integration with isotope data along trajectories

### 4. Visualization Options
Set these options in the namelist file:
```fortran
plot_bulktraj_with_humidity = True
plot_bulktraj_with_pressure = True
plot_bulktraj_with_moisture_flux = True
plot_bulktraj_with_moisturetake = True
plot_bulktraj_with_Precipitable_Water_Delta18O = True
plot_bulktraj_with_Precipitable_Water_DeltaD = True
```

## Output Files

### 1. Trajectory Files
- HYSPLIT format trajectory files
- Clipped versions for clustering analysis
- Reverse trajectories (if enabled)

### 2. NetCDF Files
- Trajectory data with meteorological parameters
- Moisture diagnostics
- Isotope compositions

### 3. Visualization Plots
- Trajectory maps with various parameters
- Moisture source regions
- Isotopic composition distributions

## Troubleshooting

Common issues and solutions:

1. HYSPLIT executable not found:
   ```bash
   # Check your namelist file hysplit path
   hysplit = ../exec
   ```

2. Missing meteorological data:
   - Ensure GDAS1 data is present in meteo_dir
   - Check file naming convention matches expected format

3. Isotope data issues:
   - Verify isotope file format matches expected structure
   - Check coordinate systems match between trajectory and isotope data

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If using this code in your research, please cite:
1. Stein, et al. (2015) NOAA's HYSPLIT Atmospheric Transport and Dispersion Modeling System, BAMS
2. The source projects (splitr, HyTraj, and PySPLIT)
3. This package

## Contact

For questions and support:
1. Open an issue on the GitHub repository
2. Check existing documentation and issues
3. Review the example namelist files in the namelist directory