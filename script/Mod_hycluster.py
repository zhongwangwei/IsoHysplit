import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pywt
from sklearn.cluster import KMeans
import cartopy.crs as ccrs
import cartopy.feature as cfeature

class HyData:
    def __init__(self, files):
        self.files = files
        # print('Reading Data....')

    def read(self):
        ds = xr.Dataset()
        ds= self._read().to_array(dim="geo")
        return ds

    def _read(self):
        import re
        files = self.files
        print("Example filename:", files[0])
        
        dates = []
        for filename in files:
            try:
                date_match = re.search(r'\d{10}', filename)
                date = pd.to_datetime(date_match.group(), format="%Y%m%d%H")
                dates.append(date)
            except:
                print(f"Error parsing date from {filename}")
                continue

        lat = pd.DataFrame(columns=dates)
        lon = pd.DataFrame(columns=dates)
        alt = pd.DataFrame(columns=dates)
        pre = pd.DataFrame(columns=dates)
        dat = [lat, lon, alt, pre]
        for filename, date in zip(files, dates):
            data = self.read_hysplit_file(filename)
            for var, col in zip(dat, data.columns):
                var[date] = data[col]
        data = xr.Dataset({"lat": lat, "lon": lon, "alt": alt, "pre": pre})
        data = data.rename({"dim_0": "step", "dim_1": "time"}).astype("float")
        return data

    def read_hysplit_file(self, filename, skp=0):
        columns = ["lat", "lon", "height", "pressure"]
        with open(filename, "r") as f:
            data = f.readlines()
        
        for num, line in enumerate(data):
            if "PRESSURE" in line:
                skp = num + 1  # Skip the header line itself
                break
        
        trajectory_data = []
        for line in data[skp + 1:]:  # Start from the line after PRESSURE
            if line.strip():
                try:
                    values = line.strip().split()
                    # The lat, lon, height, pressure values are at indices 9,10,11,12
                    trajectory_point = [float(values[9]), float(values[10]), 
                                      float(values[11]), float(values[12])]
                    trajectory_data.append(trajectory_point)
                except (IndexError, ValueError):
                    continue
        
        data = pd.DataFrame(trajectory_data, columns=columns)
        return data

class HyCluster:
    def __init__(
        self,
        data,
        scale=False,
    ):
        self.data = data
        self.scale = scale
        self.feat = HyWave(data).fit(scale=scale)

    def fit(self, kmax=50, method="KMeans", scale=False):
        labels = Trajclustering(self.feat).fit(kmax=kmax)
        self.labels = pd.DataFrame(labels).T
        return self.labels

    def get_kmeans_cluster(self, n_clusters=4):
        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(self.feat)
        labels = pd.Series(kmeans.labels_, index=self.feat.index)
        self.labels = pd.DataFrame(labels).T
        return self.labels

class HyWave:
    def __init__(
        self, data
    ):
        self.data = data
        self.time = data.time.to_pandas()

    def fit(self, scale=True):
        ln, lt = (
            self.data.sel(geo="lon").values, self.data.sel(geo="lat").values
        )
        ff = pd.concat([self._wavelet_features(lt), self._wavelet_features(ln)])
        ff.index = [
            "latmin",
            "lat25",
            "lat50",
            "lat75",
            "latmax",
            "lonmin",
            "lon25",
            "lon50",
            "lon75",
            "lonmax",
        ]
        if scale:
            ff = (ff - ff.min()) / (ff.max() - ff.min())
        return ff.T

    def _wavelet_features(self, data):
        wv = pywt.dwt(data.T, "haar")[0]
        wv = pd.DataFrame(wv, self.time).T.describe().iloc[3:]
        return wv

class Trajclustering:
    def __init__(self, data):
        self.traj = data

    def fit(self, kmax=50):
        n, wce, labels = self.get_kmeans_cluster(kmax, plot=False)
        return labels

    def _elbow_method(self, kmax=50):
        wce = []
        nums = np.arange(1, kmax)
        for num in nums:
            kmeans = KMeans(n_clusters=num, random_state=0).fit(self.traj)
            wce.append(kmeans.inertia_)

        x0, y0 = 0.0, wce[0]
        x1, y1 = float(len(wce)), wce[-1]
        elbows = []
        for index_elbow in range(1, len(wce) - 1):
            x, y = float(index_elbow), wce[index_elbow]
            segment = abs((y0 - y1) * x + (x1 - x0) * y + (x0 * y1 - x1 * y0))
            norm = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
            distance = segment / norm
            elbows.append(distance)
        n = nums[np.argmax(elbows) + 1]
        return n, wce



    def get_kmeans_cluster(self, kmax=50, plot=True):
        n, wce = self._elbow_method(kmax=kmax)
        kmeans = KMeans(n_clusters=n, random_state=0).fit(self.traj)
        labels = pd.Series(kmeans.labels_, index=self.traj.index)
        self.optim_k = n
        self.wce = wce
        if plot:
            self._plot_elbow_score(n, wce)
        return n, wce, labels

    def _plot_elbow_score(self, n, wce):
        nums = np.arange(1, len(wce) + 1)
        fig, ax = plt.subplots(1, 1, figsize=(14, 5))
        ax.plot(nums, wce, color="m")
        ax.scatter(n, wce[n - 1], color="red", marker=".", s=200)
        ax.axvline(n, ls="-.", color="k")
        ax.minorticks_on()
        ax.set_xlabel("Number of clusters")
        ax.set_ylabel("Within cluster Error")
        ax.set_title("Optimal number of clusters = %s" % n)
        plt.show()


class ClusterPlot:
    def __init__(
        self,
        data,
        cluster
    ):
        self.lat = data.sel(geo="lat").to_pandas()
        self.lon = data.sel(geo="lon").to_pandas()
        self.n_traj = len(self.lat.columns)
        self.slat = self.lat.iloc[0, 0]
        self.slon = self.lon.iloc[0, 0]
        self.cluster = cluster
        self.nclus = cluster.T.nunique().values[0]

    def mean_trajectory(self, latitudes, longitudes):
        """ Get the centroid of parcels at each timestep. """
        x = np.cos(np.radians(latitudes)) * np.cos(np.radians(longitudes))
        y = np.cos(np.radians(latitudes)) * np.sin(np.radians(longitudes))
        z = np.sin(np.radians(latitudes))

        # Get average x, y, z values
        mean_x = np.mean(x, axis=1)
        mean_y = np.mean(y, axis=1)
        mean_z = np.mean(z, axis=1)

        # Convert average values to trajectory latitudes and longitudes
        mean_longitudes = np.degrees(np.arctan2(mean_y, mean_x))
        hypotenuse = np.sqrt(mean_x ** 2 + mean_y ** 2)
        mean_latitudes = np.degrees(np.arctan2(mean_z, hypotenuse))
        return mean_latitudes, mean_longitudes

    def get_representative_trajectories(self):
        labels = self.cluster.values[0]
        columns = self.cluster.columns
        kcount = []
        clusters = ["CLUS_" + str(i + 1) for i in np.arange(self.nclus)]
        self.rep_traj_lat = pd.DataFrame(columns=clusters)
        self.rep_traj_lon = pd.DataFrame(columns=clusters)
        for num, cluster in enumerate(clusters):
            col = columns[labels == num]
            kcount.append(len(col))
            lats, lons = self.mean_trajectory(self.lat[col], self.lon[col])
            self.rep_traj_lat[cluster], self.rep_traj_lon[cluster] = (lats, lons)
        return self.rep_traj_lat, self.rep_traj_lon, kcount

    def plot_representative_trajectories(self, ax=None, cmap=plt.cm.jet, lw=3, s=200):
        xx, yy = (self.slon, self.slat)
        ax.scatter(xx, yy, color="r", s=s, transform=ccrs.PlateCarree())

        lat1, lon1, kcount = self.get_representative_trajectories()
        colors = [cmap(i) for i in np.linspace(0, 1, self.nclus)]
        for count, tr in enumerate(lat1.columns):
            lwd = lw * kcount[count] / np.sum(kcount)
            xx, yy = (lon1[tr].values, lat1[tr].values)
            ax.plot(xx, yy, color=colors[count], lw=lwd, transform=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, color="lightgrey", alpha=0.5)
        ax.add_feature(cfeature.OCEAN, color="skyblue", alpha=0.4)
        return ax


class get_cluster_Kmeans:
    def __init__(self, files, mapcorners):
        self.files = files
        self.data = HyData(files).read()
        self.labels = HyCluster(self.data).fit(kmax=10, method="KMeans")
        # Plot the clusters
        plot_corners = [mapcorners[0], mapcorners[2], mapcorners[1], mapcorners[3]]
        self.plot_clusters(plot_corners)
         
    def plot_clusters(self,mapcorners):
        # Create the map projection
        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        # Set map extent (adjust these values based on your trajectory data)
        ax.set_extent(mapcorners, crs=ccrs.PlateCarree())
        
        # Add gridlines
        ax.gridlines(draw_labels=True)
        
        # Create ClusterPlot instance and plot trajectories
        cluster_plot = ClusterPlot(self.data, self.labels)
        ax = cluster_plot.plot_representative_trajectories(ax=ax)
        
        plt.title('Trajectory Clusters')
        plt.show()
        
