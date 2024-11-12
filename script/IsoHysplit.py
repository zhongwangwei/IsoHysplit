import os
import sys
import subprocess
import re
from Mod_nml_read import NamelistReader
from Mod_traj_gen import generate_bulktraj
from Mod_traj_group import make_trajectorygroup
from Mod_traj_plot import *

def preprocess_nml(namelist_file):
    """Main function to process namelist file"""

    try:
        reader = NamelistReader()  # Create an instance first
        nml=reader.get_namelist_attr(namelist_file)
        print('namelist attributes read')
        print(nml)
        reader.write_setup_cfg(nml)
        print("Successfully created SETUP.CFG")
    except Exception as e:
        print(f"Error processing namelist file: {e}")
        sys.exit(1)
    # Convert potential list values to strings
    working_dir = nml['general']['working_dir']
    storage_dir = nml['general']['storage_dir']
    meteo_dir = nml['general']['meteo_dir']
    # If any of these are lists, take the first element
    if isinstance(working_dir, list): working_dir = working_dir[0]
    if isinstance(storage_dir, list): storage_dir = storage_dir[0]
    if isinstance(meteo_dir, list): meteo_dir = meteo_dir[0]
    
    # Convert location coordinates to float
    location = nml['general']['location']
    if isinstance(location, list):
        location = [float(coord) for coord in location]
    else:
        location = [float(coord) for coord in location.split(',')]  # If location is a comma-separated string
    #convert location to tuple
    location = tuple(location)

    #convert runtime to int
    runtime = nml['general']['runtime']
    if isinstance(runtime, list):
        runtime = int(runtime[0])
    else:
        runtime = int(runtime)
    # Convert basename to string if it's a list
    basename = nml['general']['basename']
    if isinstance(basename, list):
        basename = basename[0]    
    altitudes = nml['general']['altitudes']
    if isinstance(altitudes, list):
        altitudes = [int(alt) for alt in altitudes]  # Convert each altitude to float
    else:
        # If it's a single string, split by comma in case multiple altitudes are provided as comma-separated string
        altitudes = [int(alt.strip()) for alt in str(altitudes).split(',')]
    # Create storage directory if it doesn't exist
    os.makedirs(storage_dir, exist_ok=True)     
    nml['general']['altitudes'] = altitudes
    nml['general']['location'] = location
    nml['general']['runtime'] = runtime
    nml['general']['basename'] = basename
    nml['general']['working_dir'] = working_dir
    nml['general']['storage_dir'] = storage_dir
    nml['general']['meteo_dir'] = meteo_dir

    return nml

def main():
    if len(sys.argv) != 2:
        print("Usage: python IsoHysplit.py <namelist_file>")
        sys.exit(1)
    
    namelist_file = sys.argv[1]
    if not os.path.exists(namelist_file):
        print(f"Error: File {namelist_file} not found")
        sys.exit(1)
        
    nml = preprocess_nml(namelist_file)
    lonx = nml['general']['location'][1]
    latx = nml['general']['location'][0]
    print(f"latx: {latx}, lonx: {lonx}")
    if nml['general']['get_bulktraj']:
        print('generating bulk trajectories')
        generate_bulktraj(nml['general']['basename'], nml['general']['working_dir'], nml['general']['storage_dir'], nml['general']['meteo_dir'],
                         nml['general']['years'], nml['general']['months'], 
                         nml['general']['hours'], nml['general']['altitudes'], 
                         nml['general']['location'],  
                         nml['general']['runtime'],
                         monthslice=slice(0, 32, 1), get_reverse=True,
                         get_clipped=True, hysplit=nml['general']['hysplit'])
    if nml['general']['plot_bulktraj_with_humidity']:
        print(f"{nml['general']['storage_dir']}/{nml['general']['basename']}*")
        trajgroup = make_trajectorygroup(f"{nml['general']['storage_dir']}/{nml['general']['basename']}*")
        plot_bulktraj_with_humidity(trajgroup,nml['plot']['mapcorners'],latx=latx,lonx=lonx)
        print('plotting bulk trajectories with humidity')
    if nml['general']['plot_bulktraj_with_moisture_flux']:
        print('plotting bulk trajectories with moisture flux')
        trajgroup = make_trajectorygroup(f"{nml['general']['storage_dir']}/{nml['general']['basename']}*")
        plot_bulktraj_with_moisture_flux(trajgroup,nml['plot']['mapcorners'],latx=latx,lonx=lonx)
    if nml['general']['plot_bulktraj_with_moisturetake']:
        print('plotting bulk trajectories with moisture take')
        trajgroup = make_trajectorygroup(f"{nml['general']['storage_dir']}/{nml['general']['basename']}*")
        plot_bulktraj_with_moisturetake_new(trajgroup,nml['plot']['mapcorners'],latx=latx,lonx=lonx)
    if nml['general']['plot_bulktraj_with_Delta_D']:
        print('plotting bulk trajectories with Delta D')
        trajgroup = make_trajectorygroup(f"{nml['general']['storage_dir']}/{nml['general']['basename']}*")
        with xr.open_dataset('g2005iso-n_pwatclm.nc') as ds:
            dd=ds['dd']
        plot_bulktraj_with_Delta_D(trajgroup,dd,nml['plot']['mapcorners'],latx=latx,lonx=lonx)
if __name__ == "__main__":
    main()
