import os
import sys
import argparse
from Mod_nml_read import NamelistReader
from Mod_traj_gen import generate_bulktraj
from Mod_traj_group import make_trajectorygroup
from Mod_traj_plot import *
import pandas as pd
import xarray as xr
import glob
from Mod_hycluster import get_cluster_Kmeans

def get_or_compute_isotope_data(isotope_type, data_variable_name, file_pattern, calculation_func, isotope_dir, min_year, max_year):
    """
    Loads isotope data from a cached NetCDF file or computes it from source files.
    """
    cache_file = f"{isotope_dir}/{isotope_type}_{min_year}_{max_year}.nc"
    try:
        with xr.open_dataset(cache_file) as ds:
            print(f"Loading cached isotope data from {cache_file}")
            data = ds[data_variable_name]
    except (FileNotFoundError, OSError):
        print(f"Cached file not found or error reading it. Computing from source: {cache_file}")
        isotope_files = [
            f"{isotope_dir}/{file_pattern.format(year=year)}"
            for year in range(min_year, max_year + 1)
        ]
        print(f"Loading source isotope files: {isotope_files}")
        isotope_files.sort()

        if len(isotope_files) == 1:
            with xr.open_dataset(isotope_files[0]) as ds:
                data = calculation_func(ds).compute()
        else:
            with xr.open_mfdataset(isotope_files, parallel=True, combine="by_coords") as ds:
                data = calculation_func(ds).compute()

        data.name = data_variable_name
        data.to_netcdf(cache_file, compute=True)
        print(f"Saved computed isotope data to {cache_file}")

    return data

def main():
    parser = argparse.ArgumentParser(description="Run IsoHysplit with a given namelist file.")
    parser.add_argument('namelist_file', nargs='?', default='../namelist/main.nml', help='Path to the namelist file.')
    args = parser.parse_args()

    # Construct the absolute path to the namelist file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    namelist_path = os.path.join(script_dir, args.namelist_file)

    project_root = os.path.dirname(script_dir)
    reader = NamelistReader(project_root)
    nml = reader.get_namelist_attr(namelist_path)
    
    lonx = nml['general']['location'][1]
    latx = nml['general']['location'][0]
    print(f"latx: {latx}, lonx: {lonx}")

    if nml['general']['get_bulktraj']:
        print('generating bulk trajectories')
        os_name = os.name
        platform_system = sys.platform
        if os_name == 'nt':
            hysplit_exec = 'hyts_std.exe'
        elif platform_system == 'darwin':
            hysplit_exec = 'hyts_std_Mac'
        elif platform_system.startswith('linux'):
            hysplit_exec = 'hyts_std_Linux_x86'
        else:
            hysplit_exec = None

        nml['general']['hours'] = [f"{int(hour):02}" for hour in nml['general']['hours']]
        
        hysplit_dir = nml['general']['hysplit']
        if isinstance(hysplit_dir, list):
            hysplit_dir = hysplit_dir[0]
        hysplit_path = os.path.join(nml['general']['hysplit'], hysplit_exec)
        generate_bulktraj(
            nml['general']['basename'],
            nml['general']['working_dir'],
            nml['general']['storage_dir'],
            nml['general']['meteo_dir'],
            nml['general']['location'],
            nml['general']['runtime'],
            nml['general']['hours'],
            nml['general']['altitudes'],
            nml['general']['startdate'],
            nml['general']['enddate'],
            get_reverse=True,
            get_clipped=True,
            hysplit=hysplit_path
        )

    trajgroup = make_trajectorygroup(f"{nml['general']['storage_dir']}/{nml['general']['basename']}/traj/{nml['general']['basename']}*")
    
    min_time = min(t.data.DateTime.min() for t in trajgroup)
    max_time = max(t.data.DateTime.max() for t in trajgroup)
    min_year, max_year = min_time.year, max_time.year
    print(f"Trajectory date range: {nml['general']['startdate'][0]} to {nml['general']['enddate'][0]}")

    outdir = f"{nml['general']['storage_dir']}/{nml['general']['basename']}/"
    os.makedirs(outdir, exist_ok=True)
    
    plot_configs = {
        'plot_bulktraj_with_humidity': (plot_bulktraj_with_humidity, []),
        'plot_bulktraj_with_pressure': (plot_bulktraj_with_pressure, []),
        'plot_bulktraj_with_temperature': (plot_bulktraj_with_temperature, []),
        'plot_bulktraj_with_moisture_flux': (plot_bulktraj_with_moisture_flux, []),
        'plot_bulktraj_with_moistureuptake': (plot_bulktraj_with_moisturetake_new, []),
        'plot_bulktraj_frequency': (plot_bulktraj_frequency, [])
    }

    for key, (plot_func, args) in plot_configs.items():
        if nml['general'].get(key):
            plot_func(trajgroup, outdir, nml['plot']['mapcorners'], latx=latx, lonx=lonx, *args)

    isotope_dir = nml['general']['isotope_dir']
    map_corners = nml['plot']['mapcorners']

    isotope_tasks = {
        'Precipitable_Water_Delta18O': ('do', 'g{year}iso-n.nc', lambda ds: (ds['pwat1clm'] / ds['pwatclm'] - 1) * 1000, plot_bulktraj_with_Delta, ('Precipitable_Water', '18O')),
        'Precipitable_Water_DeltaD': ('dd', 'g{year}iso-n.nc', lambda ds: (ds['pwat2clm'] / ds['pwatclm'] - 1) * 1000, plot_bulktraj_with_Delta, ('Precipitable_Water', 'D')),
        'Precipitation_Water_Delta18O': ('do', 'g{year}iso-n.nc', lambda ds: (ds['prate1'] / ds['prate'] - 1) * 1000, plot_bulktraj_with_Delta, ('Precipitation_Water', '18O')),
        'Precipitation_Water_DeltaD': ('dd', 'g{year}iso-n.nc', lambda ds: (ds['prate2'] / ds['prate'] - 1) * 1000, plot_bulktraj_with_Delta, ('Precipitation_Water', 'D')),
        'Vapor_Water_Delta18O': ('do', 'g{year}iso-n_pgb.nc', lambda ds: (ds['spfh1prs'] / ds['spfhprs'] - 1) * 1000, plot_bulktraj_with_Delta_with_level, ('Vapor_Water', '18O')),
        'Vapor_Water_DeltaD': ('dd', 'g{year}iso-n_pgb.nc', lambda ds: (ds['spfh2prs'] / ds['spfhprs'] - 1) * 1000, plot_bulktraj_with_Delta_with_level, ('Vapor_Water', 'D')),
        'Vapor_Water_D_excess': ('dex', 'g{year}iso-n_pgb.nc', lambda ds: ((ds['spfh2prs'] / ds['spfhprs'] - 1) * 1000) - 8 * ((ds['spfh1prs'] / ds['spfhprs'] - 1) * 1000), plot_bulktraj_with_Delta_with_level, ('Vapor_Water', 'D_excess')),
        'Precipitation_Water_D_excess': ('dex', 'g{year}iso-n.nc', lambda ds: ((ds['prate2'] / ds['prate'] - 1) * 1000) - 8 * ((ds['prate1'] / ds['prate'] - 1) * 1000), plot_bulktraj_with_Delta, ('Precipitation_Water', 'D_excess')),
        'Precipitable_Water_D_excess': ('dex', 'g{year}iso-n.nc', lambda ds: ((ds['pwat2clm'] / ds['pwatclm'] - 1) * 1000) - 8 * ((ds['pwat1clm'] / ds['pwatclm'] - 1) * 1000), plot_bulktraj_with_Delta, ('Precipitable_Water', 'D_excess')),
    }

    for plot_key, (var_name, file_pattern, calc_func, plot_func, plot_args) in isotope_tasks.items():
        if nml['general'].get(f'plot_bulktraj_with_{plot_key}'):
            data = get_or_compute_isotope_data(plot_key, var_name, file_pattern, calc_func, isotope_dir, min_year, max_year)
            plot_func(trajgroup, data, *plot_args, outdir=outdir, mapcorners=map_corners, latx=latx, lonx=lonx)

    if nml['general'].get('get_cluster_Kmeans'):
        files = sorted(glob.glob(f"{nml['general']['storage_dir']}/{nml['general']['basename']}/traj/{nml['general']['basename']}*"))
        get_cluster_Kmeans(files, mapcorners=map_corners, save_figure=True, figure_path='trajectory_clusters.png', save_data=True, data_path=outdir)

if __name__ == "__main__":
    main()