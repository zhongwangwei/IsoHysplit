import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from Mod_traj_plotlib import *
import os
from joblib import Parallel, delayed

def plot_bulktraj_with_humidity(trajgroup,outdir,mapcorners,latx=1.352,lonx=103.820):
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
   bmap_params = MapDesign(mapcorners, standard_pm,latx=latx,lonx=lonx)

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
   output_dir = f'{outdir}/plot_bulktraj_with_humidity'
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
   ##plt.show()

def plot_bulktraj_with_pressure(trajgroup,outdir,mapcorners,latx=1.352,lonx=103.820):
   fig,ax0 = plt.subplots(nrows=1, figsize=(10,7))
   standard_pm = None
   colors = np.linspace(0.10, 0.9, trajgroup.trajcount)
   bmap_params = MapDesign(mapcorners, standard_pm,latx=latx,lonx=lonx)
   map_scatter  = bmap_params.make_basemap(ax=ax0)
   allpressure = []
   for traj in trajgroup:
       allpressure.extend(traj.data.Pressure.astype(np.float64).values)
    
   vmin = np.nanpercentile(allpressure, 5)
   vmax = np.nanpercentile(allpressure, 95)
   # Round vmin and vmax to nearest 10
   vmin = np.floor(vmin / 10) * 10
   vmax = np.ceil(vmax / 10) * 10

   for traj in trajgroup[::1]:
      mappable = traj_scatter(
         traj.data.Pressure.astype(np.float64).values,
         traj.data.geometry.apply(lambda p: p.x).values,
         traj.data.geometry.apply(lambda p: p.y).values,
         map_scatter, colormap=plt.cm.viridis, # cnormalize='sqrt',
         vmin=vmin, vmax=vmax, size=3,suppress_printmsg=True)

   # Make colorbar on its own axis
   cax_position = [0.2, 0.1, 0.6, 0.05]
   cax, cbar = make_cax_cbar(fig, cax_position, mappable,
   tick_fs=14, label_fs=16, cbar_label='Pressure (hPa)',
                                    labelpad=12);
   ticklabels = cbar.ax.get_xticklabels()
   newlabels = []
   # Create output directory if it doesn't exist
   output_dir = f'{outdir}/plot_bulktraj_with_pressure'
   os.makedirs(output_dir, exist_ok=True)
   # Save the plot
   plt.savefig(f'{output_dir}/plot_bulktraj_with_pressure.png', bbox_inches='tight', dpi=300)


def plot_bulktraj_with_moisture_flux(trajgroup,outdir,mapcorners,latx=1.352,lonx=103.820):
   for traj in trajgroup:
      traj.calculate_distance()
      traj.calculate_vector()
      traj.calculate_moistureflux()

   fig,ax0 = plt.subplots(nrows=1, figsize=(10,7))
   standard_pm = None

   """
   Plotting ``Trajectory`` Paths
   -----------------------------
   For this example, we will color-code by initialization (t=0) altitude,
   (500, 1000, or 1500 m), which can be accessed via ``Trajectory.data.geometry``,
   a ``GeoSeries`` of Shapely ``Point`` objects.

   We can store the trajectory color in ``Trajectory.trajcolor`` for convenience.

   """
   colors = np.linspace(0.10, 0.9, trajgroup.trajcount)

   bmap_params = MapDesign(mapcorners, standard_pm,latx=latx,lonx=lonx)

   map_scatter  = bmap_params.make_basemap(ax=ax0)

   min_mf = 0
   max_mf = 200
   for traj in trajgroup:
      if traj.data.Moisture_Flux.max() > max_mf:
         max_mf = traj.data.Moisture_Flux.max()
      if traj.data.Moisture_Flux.min() < min_mf:
         min_mf = traj.data.Moisture_Flux.min()

   for traj in trajgroup[::1]:
      mappable = traj_scatter(
         traj.data.Moisture_Flux.astype(np.float64).values,
         traj.data.geometry.apply(lambda p: p.x).values,
         traj.data.geometry.apply(lambda p: p.y).values,
         map_scatter, colormap=plt.cm.viridis, # cnormalize='sqrt',
         vmin=min_mf, vmax=200.0, size=3,suppress_printmsg=True)
   # Make colorbar on its own axis
   cax_position = [0.2, 0.1, 0.6, 0.05]
   cax, cbar = make_cax_cbar(fig, cax_position, mappable,
      tick_fs=14, label_fs=16, cbar_label='Moisture Flux ((g/kg)(m/s))',
                                 labelpad=12)
   ticklabels = cbar.ax.get_xticklabels()
   newlabels = []
   output_dir = f'{outdir}/Moisture_Flux'
   os.makedirs(output_dir, exist_ok=True)

   plt.savefig(f'{output_dir}/Moisture_Flux.png', bbox_inches='tight', dpi=300)

def plot_bulktraj_with_moisturetake(trajgroup,outdir, mapcorners,latx=1.352,lonx=103.820):
   import matplotlib.pyplot as plt
   import numpy as np
   import xarray as xr
   import pandas as pd
   from matplotlib import colors
   import cartopy .crs as ccrs
   import cartopy.feature as cfeature
   from matplotlib import cm
   from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
   import cartopy.feature as cfeature
   import matplotlib
   from pylab import rcParams

   def find_nearest(array, value, k):
      array = np.asarray(array)
      idx = np.argsort(abs(array - value))[:k]
      return idx


   Intime=pd.date_range("2000-01-15",freq="1M",periods=1)
   lat=np.arange(-89.5, 90.5, 1.0)
   lon=np.arange(-179.5, 180.5, 1.0)
   dq =np.zeros((1,180,360))
   ds = xr.Dataset({'dq': (('time','lat','lon'), dq),},
                  coords={'lon': lon, 'lat':lat,'time':(('time'),Intime)})

   moisture = []
   for traj in trajgroup:
      traj.calculate_distance()
      traj.calculate_vector()
      traj.calculate_moistureflux()
      traj.load_reversetraj()
      traj.calculate_integrationerr()
      traj.moisture_uptake(precipitation=-0.2,evaporation=0.2,interval=6)
      moisture.append(traj)
   #print(traj.uptake.below.astype(np.float64).values)

   min_above_total = 0.0
   max_above_total = 1.0

   for traj in moisture:
      if traj.uptake.above_total.max() > max_above_total:
            max_above_total = traj.uptake.dq.max()
            #print(max_above_total)
      if traj.uptake.above_total.min() < min_above_total:
            min_above_total = traj.uptake.above_total.min()
            #print(min_above_total)
   for traj in moisture:
      x=traj.uptake.geometry.apply(lambda p: p.x).values
      y=traj.uptake.geometry.apply(lambda p: p.y).values
      for i, ix ,iy in zip(range(len(x)),x,y):
            iix=find_nearest(lon, ix, 1)
            iiy=find_nearest(lat, iy, 1)
            ds.dq[0,iiy,iix]=ds.dq[0,iiy,iix]+np.nan_to_num(traj.uptake.dq.astype(np.float64).values[i].flatten())
   ds1=ds.dq/ds.dq.sum()
   out_dir = f'{outdir}/Moisture_Take'
   os.makedirs(out_dir, exist_ok=True)
   ds.to_netcdf(f'{out_dir}/moisturetake.nc',engine='netcdf4')
   ds1.to_netcdf(f'{out_dir}/moisturetake_fraction.nc',engine='netcdf4')
   ds=ds.where(ds>0.01,drop=True)
   fig = plt.figure()
   ef, ax = plt.subplots(1,1,figsize=(10,5),subplot_kw={'projection': ccrs.PlateCarree()})

   # Convert mapcorners to float if they are strings
   mapcorners = [float(x) for x in mapcorners]

   ax.set_extent(mapcorners, crs=ccrs.PlateCarree())
   ax.add_feature(cfeature.LAND)
   ax.add_feature(cfeature.OCEAN)
   ax.add_feature(cfeature.COASTLINE)
   #ax.add_feature(cfeature.BORDERS, linestyle=':')
   ax.add_feature(cfeature.LAKES, alpha=0.5)
   ax.add_feature(cfeature.RIVERS)

   # Generate ticks based on mapcorners [lonmin, lonmax, latmin, latmax]
   lon_ticks = np.arange(np.floor(mapcorners[0]/10)*10, 
                        np.ceil(mapcorners[1]/10)*10, 10)
   lat_ticks = np.arange(np.floor(mapcorners[2]/5)*5,
                        np.ceil(mapcorners[3]/5)*5, 5)
   
   ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
   ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())
   lon_formatter = LongitudeFormatter(zero_direction_label=True)
   lat_formatter = LatitudeFormatter()
   ax.xaxis.set_major_formatter(lon_formatter)
   ax.yaxis.set_major_formatter(lat_formatter)

   levels=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
   levels=[0, 2, 4,   6,  8, 10, 12, 14, 16, 18, 20, 22, 24]
   #ds=xr.open_dataset("test.nc",decode_times=False)
   ds.dq.plot(ax=ax, levels=levels, cbar_kwargs={"label": "Accumulated specific humidity (g/kg)"})

   ax.set_title("")
   ax.set_xlabel("")
   # Make it nice
   plt.tight_layout()
   # Save the plot
   plt.savefig(f'{out_dir}/moisturetake.png', bbox_inches='tight', dpi=300)
   del(ds,ds1)

def plot_bulktraj_with_startpoint_endpoint(trajgroup,outdir,mapcorners,latx=1.352,lonx=103.820):
   import numpy as np
   import xarray as xr
   import pandas as pd
   def find_nearest(array, value, k):
      array = np.asarray(array)
      idx = np.argsort(abs(array - value))[:k]
      return idx
   def traj_moisture_uptake(trajgroup):
      lat=np.arange(-89.875, 90, 0.25)
      lon=np.arange(-179.875, 180, 0.25)
      dq =np.zeros((1,720,1440))
      #Intime=datetime.datetime(year=int(year),month=int(month),day=15)
      Intime=pd.date_range("2000-01-15",freq="1M",periods=1)
      lat=np.arange(-89.5, 90.5, 1.0)
      lon=np.arange(-179.5, 180.5, 1.0)
      sq =np.zeros((1,180,360))
      dq =np.zeros((1,180,360))

      ds = xr.Dataset({'sq': (('time','lat','lon'), sq),},
                  coords={'lon': lon, 'lat':lat,'time':(('time'),Intime)})

      moisture = []
      iii=1
      k2=0.0
      k3=0.0
      kk=0.0
      for traj in trajgroup:
         #try:
         x=traj.data.geometry.apply(lambda p: p.x).values[-1]
         y=traj.data.geometry.apply(lambda p: p.y).values[-1]
         k=traj.data.Specific_Humidity.astype(np.float64).values[-1]
         kk=traj.data.Specific_Humidity.astype(np.float64).values[0]
         iix=find_nearest(lon, x, 1)
         iiy=find_nearest(lat, y, 1)
         ds.sq[0,iiy,iix]=ds.sq[0,iiy,iix]+k
         k2=k2+k
         k3=k3+kk
         iii=iii+1
      ss=k2/iii
      ss3=k3/iii
      for traj in trajgroup:
         traj.calculate_distance()
         traj.calculate_vector()
         traj.calculate_moistureflux()
         traj.load_reversetraj()
         traj.calculate_integrationerr()
         traj.moisture_uptake(precipitation=-0.05,evaporation=0.05,interval=3)
         moisture.append(traj)
      k3=0.0
      for traj in moisture:
         x=traj.uptake.geometry.apply(lambda p: p.x).values
         y=traj.uptake.geometry.apply(lambda p: p.y).values
         for i, ix ,iy in zip(range(len(x)),x,y):
            iix=find_nearest(lon, ix, 1)
            iiy=find_nearest(lat, iy, 1)
            dq[0,iiy,iix]=dq[0,iiy,iix]+np.nan_to_num(traj.uptake.dq.astype(np.float64).values[i].flatten())
            k3=np.nan_to_num(traj.uptake.dq.astype(np.float64).values[i].flatten())+k3
      ss2=k3/iii
      # Create dataset with moisture uptake, source and endpoint values
      ds_moisture = xr.Dataset(
          {
              'uptake': (('time', 'lat', 'lon'), dq),
              'source': ss2,
              'endpoint': ss3
          },
          coords={
              'lon': lon,
              'lat': lat,
              'time': (('time'), Intime)
          }
      )
      
      # Save to netCDF file
      out_dir = f'{outdir}/Moisture_startpoint_endpoint'
      os.makedirs(out_dir, exist_ok=True)
      ds_moisture.to_netcdf(f'{out_dir}/Moisture_startpoint_endpoint.nc')
      return None

def plot_bulktraj_with_moisturetake_new(trajgroup,outdir,mapcorners,latx=1.352,lonx=103.820):
   import matplotlib.pyplot as plt
   import numpy as np
   import xarray as xr
   import pandas as pd
   from matplotlib import colors
   import cartopy .crs as ccrs
   import cartopy.feature as cfeature
   from matplotlib import cm
   from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
   import cartopy.feature as cfeature
   import matplotlib
   from pylab import rcParams

   def find_nearest(array, value, k):
      array = np.asarray(array)
      idx = np.argsort(abs(array - value))[:k]
      return idx


   Intime=pd.date_range("2000-01-15",freq="1M",periods=1)
   lat=np.arange(-89, 90.5, 2.0)
   lon=np.arange(-179, 180.5, 2.0)
   dq =np.zeros((1,90,180))
   ds = xr.Dataset({'dq': (('time','lat','lon'), dq),},
                  coords={'lon': lon, 'lat':lat,'time':(('time'),Intime)})

   moisture = []
   for traj in trajgroup:
      traj.calculate_distance()
      traj.calculate_vector()
      traj.calculate_moistureflux()
      traj.load_reversetraj()
      traj.calculate_integrationerr()
      traj.moisture_uptake(precipitation=-0.2,evaporation=0.2,interval=6)
      moisture.append(traj)
   #print(traj.uptake.below.astype(np.float64).values)

   min_above_total = 0.0
   max_above_total = 1.0

   for traj in moisture:
      if traj.uptake.above_total.max() > max_above_total:
            max_above_total = traj.uptake.dq.max()
            #print(max_above_total)
      if traj.uptake.above_total.min() < min_above_total:
            min_above_total = traj.uptake.above_total.min()
            #print(min_above_total)
   for traj in moisture:
      x=traj.uptake.geometry.apply(lambda p: p.x).values
      y=traj.uptake.geometry.apply(lambda p: p.y).values
      for i, ix ,iy in zip(range(len(x)),x,y):
            iix=find_nearest(lon, ix, 1)
            iiy=find_nearest(lat, iy, 1)
            ds.dq[0,iiy,iix]=ds.dq[0,iiy,iix]+np.nan_to_num(traj.uptake.dq.astype(np.float64).values[i].flatten())
   ds1=ds.dq/ds.dq.sum()
   out_dir = f'{outdir}/Moisture_Take'
   os.makedirs(out_dir, exist_ok=True)
   ds.to_netcdf(f'{out_dir}/moisturetake.nc',engine='netcdf4')
   ds1.to_netcdf(f'{out_dir}/moisturetake_fraction.nc',engine='netcdf4')
   #ds=ds.where(ds>0.01,drop=True)
   fig, ax0 = plt.subplots(nrows=1, figsize=(10,7))
   standard_pm = None
   bmap_params = MapDesign(mapcorners, standard_pm, latx=latx, lonx=lonx)
   map_scatter = bmap_params.make_basemap(ax=ax0)
   
   # Plot the data using pcolormesh instead of contourf to avoid interpolation
   levels = [0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72]
   mappable = map_scatter.pcolormesh(ds.lon, ds.lat, ds.dq[0], 
                                    vmin=min(levels), 
                                    vmax=max(levels),
                                    cmap='viridis')
   
   # Add colorbar
   cax_position = [0.2, 0.1, 0.6, 0.05]
   cax, cbar = make_cax_cbar(fig, cax_position, mappable,
       tick_fs=14, label_fs=16, 
       cbar_label='Accumulated specific humidity (g/kg)',
       labelpad=12)
   


   plt.savefig(f'{out_dir}/moisturetake_new.png', bbox_inches='tight', dpi=300)
   del(ds,ds1)

def plot_bulktraj_with_Delta_D(trajgroup,delta,mapcorners,latx=1.352,lonx=103.820):
   def process_point(ix, iy, tim, delta):
      # Find nearest time, lon, lat values in delta dataset
      time_idx = delta.time.sel(time=tim, method="nearest")
      try:
         lat_idx = delta.lat.sel(lat=iy, method="nearest")
         lon_idx = delta.lon.sel(lon=ix, method="nearest")
         delta_d_value = delta.sel(time=time_idx, lat=lat_idx, lon=lon_idx).values
      except:
         lat_idx = delta.latitude.sel(latitude=iy, method="nearest")
         lon_idx = delta.longitude.sel(longitude=ix, method="nearest")
         delta_d_value = delta.sel(time=time_idx, latitude=lat_idx, longitude=lon_idx).values
      return delta_d_value

   # Add trajectory data
   for i, traj in enumerate(trajgroup):
      # Extract coordinates and data
      lons = traj.data.geometry.apply(lambda p: p.x).values
      lats = traj.data.geometry.apply(lambda p: p.y).values
      times = traj.data.DateTime.values
      
      # Create a list to store delta_d values
      delta_d_values = Parallel(n_jobs=-1)(
          delayed(process_point)(ix, iy, tim, delta)
          for ix, iy, tim in zip(lons, lats, times)
      )
      
      # Add the delta_d values as a new column to traj.data
      traj.data['Delta_D'] = delta_d_values
   
   # Create dataset from trajectory group
   ds = xr.Dataset()
   
   # Add trajectory data
   for i, traj in enumerate(trajgroup):
       # Extract coordinates and data
       lons = traj.data.geometry.apply(lambda p: p.x).values
       lats = traj.data.geometry.apply(lambda p: p.y).values
       delta_d = traj.data.Delta_D.astype(np.float64).values
       
       # Add to dataset with trajectory number as dimension
       ds[f'trajectory_{i}_delta_d'] = xr.DataArray(
           data=delta_d,
           dims=['Timestep'],
           coords={'Timestep': traj.data.index,
                  'latitude': ('Timestep', lats),
                  'longitude': ('Timestep', lons)}
       )

   # Create output directory if it doesn't exist
   output_dir = '../output/Delta_D'
   os.makedirs(output_dir, exist_ok=True)

   # Save to netCDF file
   ds.to_netcdf(f'{output_dir}/trajectories_with_delta_d.nc')

   # Add plotting section
   fig, ax0 = plt.subplots(nrows=1, figsize=(10,7))
   standard_pm = None
   bmap_params = MapDesign(mapcorners, standard_pm, latx=latx, lonx=lonx)
   map_scatter = bmap_params.make_basemap(ax=ax0)

   for traj in trajgroup[::1]:
      mappable = traj_scatter(
         traj.data.Delta_D.astype(np.float64).values,
         traj.data.geometry.apply(lambda p: p.x).values,
         traj.data.geometry.apply(lambda p: p.y).values,
         map_scatter, colormap=plt.cm.viridis,
         vmin=-200.0, vmax=-80.0, size=3, suppress_printmsg=True)

   # Make colorbar on its own axis
   cax_position = [0.2, 0.1, 0.6, 0.05]
   cax, cbar = make_cax_cbar(fig, cax_position, mappable,
      tick_fs=14, label_fs=16, cbar_label='δD (‰)',
      labelpad=12)

   # Save the plot
   output_dir = '../output/Delta_D'
   plt.savefig(f'{output_dir}/trajectory_delta_d.png', bbox_inches='tight', dpi=300)
   #plt.show()

def plot_bulktraj_with_Delta(trajgroup,delta,varname,isotope_type,outdir,mapcorners,latx=1.352,lonx=103.820):
   def process_point(ix, iy, tim, delta):
      # Find nearest time, lon, lat values in delta dataset
      time_idx = delta.time.sel(time=tim, method="nearest")
      try:
         lat_idx = delta.lat.sel(lat=iy, method="nearest")
         lon_idx = delta.lon.sel(lon=ix, method="nearest")
         delta_value = delta.sel(time=time_idx, lat=lat_idx, lon=lon_idx).values
      except:
         lat_idx = delta.latitude.sel(latitude=iy, method="nearest")
         lon_idx = delta.longitude.sel(longitude=ix, method="nearest")
         delta_value = delta.sel(time=time_idx, latitude=lat_idx, longitude=lon_idx).values
      return delta_value

   # Add trajectory data
   for i, traj in enumerate(trajgroup):
      # Extract coordinates and data
      lons = traj.data.geometry.apply(lambda p: p.x).values
      lats = traj.data.geometry.apply(lambda p: p.y).values
      times = traj.data.DateTime.values
      
      # Create a list to store delta_d values
      delta_values = Parallel(n_jobs=-1)(
          delayed(process_point)(ix, iy, tim, delta)
          for ix, iy, tim in zip(lons, lats, times)
      )
      
      # Add the delta_d values as a new column to traj.data
      traj.data[f'Delta_{isotope_type}'] = delta_values
   
   # Create dataset from trajectory group
   ds = xr.Dataset()
   
   # Add trajectory data
   for i, traj in enumerate(trajgroup):
       # Extract coordinates and data
       lons = traj.data.geometry.apply(lambda p: p.x).values
       lats = traj.data.geometry.apply(lambda p: p.y).values
       delta = traj.data[f'Delta_{isotope_type}'].astype(np.float64).values
       
       # Add to dataset with trajectory number as dimension
       ds[f'trajectory_{i}_delta_{isotope_type}'] = xr.DataArray(
           data=delta,
           dims=['Timestep'],
           coords={'Timestep': traj.data.index,
                  'latitude': ('Timestep', lats),
                  'longitude': ('Timestep', lons)}
       )

   # Create output directory if it doesn't exist
   output_dir = f'{outdir}/Delta_{varname}_{isotope_type}'
   os.makedirs(output_dir, exist_ok=True)

   # Save to netCDF file
   ds.to_netcdf(f'{output_dir}/trajectories_with_Delta_{varname}_{isotope_type}.nc')

   # Calculate vmin and vmax based on 5th and 95th percentiles across all trajectories
   all_delta_values = []
   for traj in trajgroup:
       all_delta_values.extend(traj.data[f'Delta_{isotope_type}'].astype(np.float64).values)
    
   vmin = np.nanpercentile(all_delta_values, 5)
   vmax = np.nanpercentile(all_delta_values, 95)
   # Round vmin and vmax to nearest 10
   vmin = np.floor(vmin / 10) * 10
   vmax = np.ceil(vmax / 10) * 10

   # Add plotting section
   fig, ax0 = plt.subplots(nrows=1, figsize=(10,7))
   standard_pm = None
   bmap_params = MapDesign(mapcorners, standard_pm, latx=latx, lonx=lonx)
   map_scatter = bmap_params.make_basemap(ax=ax0)

   #set vmin and vmax for the rest of the trajectories
   for traj in trajgroup[::1]:
      mappable = traj_scatter(
         traj.data[f'Delta_{isotope_type}'].astype(np.float64).values,
         traj.data.geometry.apply(lambda p: p.x).values,
         traj.data.geometry.apply(lambda p: p.y).values,
         map_scatter, colormap=plt.cm.viridis,
         vmin=vmin, vmax=vmax, size=3, suppress_printmsg=True)

   # Make colorbar on its own axis
   cax_position = [0.2, 0.1, 0.6, 0.05]
   cax, cbar = make_cax_cbar(fig, cax_position, mappable,
      tick_fs=14, label_fs=16, cbar_label=f'δ^{isotope_type} (‰)',
      labelpad=12)

   # Save the plot
   plt.savefig(f'{output_dir}/trajectory_delta_{varname}_{isotope_type}.png', bbox_inches='tight', dpi=300)
   #plt.show()

def plot_bulktraj_with_Delta_with_level(trajgroup,delta,varname,isotope_type,outdir,mapcorners,latx=1.352,lonx=103.820):
   def process_point(ix, iy, plev, tim, delta):
      # Find nearest time, lon, lat values in delta dataset
      time_idx = delta.time.sel(time=tim, method="nearest")
      try:
         lat_idx = delta.lat.sel(lat=iy, method="nearest")
         lon_idx = delta.lon.sel(lon=ix, method="nearest")
         level_idx = delta.levels.sel(levels=plev, method="nearest")
         delta_value = delta.sel(time=time_idx, lat=lat_idx, lon=lon_idx, levels=level_idx).values
      except:
         lat_idx = delta.latitude.sel(latitude=iy, method="nearest")
         lon_idx = delta.longitude.sel(longitude=ix, method="nearest")
         level_idx = delta.levels.sel(levels=plev, method="nearest")
         delta_value = delta.sel(time=time_idx, latitude=lat_idx, longitude=lon_idx, levels=level_idx).values
      return delta_value

   # Add trajectory data
   for i, traj in enumerate(trajgroup):
      # Extract coordinates and data
      lons = traj.data.geometry.apply(lambda p: p.x).values
      lats = traj.data.geometry.apply(lambda p: p.y).values
      times = traj.data.DateTime.values
      pressure = traj.data.Pressure.astype(np.float64).values
      # Create a list to store delta_d values
      delta_values = Parallel(n_jobs=-1)(
          delayed(process_point)(ix, iy, plev, tim, delta)
          for ix, iy, plev, tim in zip(lons, lats, pressure, times)
      )
      
      # Add the delta_d values as a new column to traj.data
      traj.data[f'Delta_{isotope_type}'] = delta_values
   
   # Create dataset from trajectory group
   ds = xr.Dataset()
   
   # Add trajectory data
   for i, traj in enumerate(trajgroup):
      # Extract coordinates and data
      lons = traj.data.geometry.apply(lambda p: p.x).values
      lats = traj.data.geometry.apply(lambda p: p.y).values
      delta = traj.data[f'Delta_{isotope_type}'].astype(np.float64).values
      
      # Add to dataset with trajectory number as dimension
      ds[f'trajectory_{i}_delta_{isotope_type}'] = xr.DataArray(
         data=delta,
         dims=['Timestep'],
         coords={'Timestep': traj.data.index,
               'pressure': ('Timestep', pressure),
               'latitude': ('Timestep', lats),
               'longitude': ('Timestep', lons)}
      )

   # Create output directory if it doesn't exist
   output_dir = f'{outdir}/Delta_{varname}_{isotope_type}'
   os.makedirs(output_dir, exist_ok=True)

   # Save to netCDF file
   ds.to_netcdf(f'{output_dir}/trajectories_with_Delta_{varname}_{isotope_type}.nc')

   # Calculate vmin and vmax based on 5th and 95th percentiles across all trajectories
   all_delta_values = []
   for traj in trajgroup:
       all_delta_values.extend(traj.data[f'Delta_{isotope_type}'].astype(np.float64).values)
    
   vmin = np.nanpercentile(all_delta_values, 5)
   vmax = np.nanpercentile(all_delta_values, 95)
   # Round vmin and vmax to nearest 10
   vmin = np.floor(vmin / 10) * 10
   vmax = np.ceil(vmax / 10) * 10

   # Add plotting section
   fig, ax0 = plt.subplots(nrows=1, figsize=(10,7))
   standard_pm = None
   bmap_params = MapDesign(mapcorners, standard_pm, latx=latx, lonx=lonx)
   map_scatter = bmap_params.make_basemap(ax=ax0)

   #set vmin and vmax for the rest of the trajectories
   for traj in trajgroup[::1]:
      mappable = traj_scatter(
         traj.data[f'Delta_{isotope_type}'].astype(np.float64).values,
         traj.data.geometry.apply(lambda p: p.x).values,
         traj.data.geometry.apply(lambda p: p.y).values,
         map_scatter, colormap=plt.cm.viridis,
         vmin=vmin, vmax=vmax, size=3, suppress_printmsg=True)

   # Make colorbar on its own axis
   cax_position = [0.2, 0.1, 0.6, 0.05]
   cax, cbar = make_cax_cbar(fig, cax_position, mappable,
      tick_fs=14, label_fs=16, cbar_label=f'δ^{isotope_type} (‰)',
      labelpad=12)

   # Save the plot
   plt.savefig(f'{output_dir}/trajectory_delta_{varname}_{isotope_type}.png', bbox_inches='tight', dpi=300)
   #plt.show()
