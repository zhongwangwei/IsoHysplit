import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from Mod_traj_plotlib import *

def plot_bulktraj_with_humidity(trajgroup,mapcorners):
   fig,ax0 = plt.subplots(nrows=1, figsize=(10,7))
   #param_dict = {'projection':'lcc', 'latlon_labelspacing':(10,30),
   #               'latlon_spacing':(10,15), 'latlon_fs':14}
   standard_pm = None

   """
   Plotting ``Trajectory`` Paths
   -----------------------------
   For this example, we will color-code by initialization (t=0) altitude,
   (500, 1000, or 1500 m), which can be accessed via ``Trajectory.data.geometry``,
   a ``GeoSeries`` of Shapely ``Point`` objects.

   We can store the trajectory color in ``Trajectory.trajcolor`` for convenience.

   """
   #color_dict = {
   #              1000.0 : 'orange'
   #              }
   colors = np.linspace(0.10, 0.9, trajgroup.trajcount)

   #bmap_params = pysplit.MapDesign(mapcorners, standard_pm,** param_dict)
   bmap_params = MapDesign(mapcorners, standard_pm)

   """
   Once the ``MapDesign`` is initialized it can be used to draw a map:

   """
   map_scatter  = bmap_params.make_basemap(ax=ax0)
   #for traj in trajgroup:
   #    altitude0 = traj.data.geometry.apply(lambda p: p.z)[0]
   #    traj.trajcolor = color_dict[altitude0]

   """
   For display purposes, let's plot only every fifth ``Trajectory``.  The lats,
   lons are obtained by unpacking the ``Trajectory.Path``
   (Shapely ``LineString``) xy coordinates.

   """

   for traj in trajgroup[::1]:
      mappable = traj_scatter(
         traj.data.Specific_Humidity.astype(np.float64).values,
         traj.data.geometry.apply(lambda p: p.x).values,
         traj.data.geometry.apply(lambda p: p.y).values,
         map_scatter, colormap=plt.cm.viridis, # cnormalize='sqrt',
         vmin=0.0, vmax=20.0, size=3,suppress_printmsg=True)
   # Make colorbar on its own axis
   cax_position = [0.2, 0.1, 0.6, 0.05]
   cax, cbar = make_cax_cbar(fig, cax_position, mappable,
   tick_fs=14, label_fs=16, cbar_label='Specific Humidity ((g/kg))',
                                    labelpad=12);
   ticklabels = cbar.ax.get_xticklabels()
   newlabels = []
   # Create output directory if it doesn't exist
   output_dir = '../output/plot_bulktraj_with_humidity'
   os.makedirs(output_dir, exist_ok=True)

   # Save the plot
   plt.savefig(f'{output_dir}/trajectory_plot.png', bbox_inches='tight', dpi=300)

   
   # Create dataset from trajectory group
   ds = xr.Dataset()
   
   # Add trajectory data
   for i, traj in enumerate(trajgroup):
       # Extract coordinates and data
       lons = traj.data.geometry.apply(lambda p: p.x).values
       lats = traj.data.geometry.apply(lambda p: p.y).values
       spec_hum = traj.data.Specific_Humidity.astype(np.float64).values
       pressure = traj.data.Pressure.astype(np.float64).values
       
       # Add to dataset with trajectory number as dimension
       ds[f'trajectory_{i}_humidity'] = xr.DataArray(
           data=spec_hum,
           dims=['Timestep'],
           coords={'Timestep': traj.data.index,
                  'latitude': ('Timestep', lats),
                  'longitude': ('Timestep', lons)}
       )
       
       ds[f'trajectory_{i}_pressure'] = xr.DataArray(
           data=pressure,
           dims=['Timestep'],
           coords={'Timestep': traj.data.index,
                  'latitude': ('Timestep', lats),
                  'longitude': ('Timestep', lons)}
       )

   # Save to netCDF file
   ds.to_netcdf(f'{output_dir}/Bulktraj_with_humidity.nc')
   plt.show()

def plot_bulktraj_with_moisture_flux(trajgroup,mapcorners):
   pass
