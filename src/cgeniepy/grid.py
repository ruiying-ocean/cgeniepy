import numpy as np
import xarray as xr
from importlib.resources import files
import pandas as pd

from scipy.interpolate import (
    RegularGridInterpolator,
    LinearNDInterpolator,
    NearestNDInterpolator,
    CubicSpline
)

class GridOperation:

    """A set of operations on grid/coordinate data
    """

    def get_genie_lon(self, N=36, edge=False, offset_start=-260):
        """get GENIE longitude in 10 degree resolution, if edge is False, then return midpoint

        :param N: number of grid points
        :param edge: if True, return edge points
        :param offset_start: the par_grid_lon_offset option in main configuration file
        """
        resolution = 360 / N
        if edge:
            lon_edge = np.linspace(offset_start, 360+offset_start, N + 1)
            return lon_edge
        else:
            lon = np.linspace(offset_start+resolution/2, 360+offset_start-resolution/2, N)
            return lon

    def get_genie_lat(self, N=36, edge=False):
        """
        return cGENIE latitude in log-sine normally degree resolution,
        if edge is False, then return midpoint
        """
        if edge:
            lat_edge = np.rad2deg(np.arcsin(np.linspace(-1, 1, N + 1)))
            return lat_edge
        else:
            lat_min = -1 + (1 / N)
            lat_max = 1 - (1 / N)
            lat = np.rad2deg(np.arcsin(np.linspace(lat_min, lat_max, N)))
            return lat

    def get_genie_depth(self, N=16, edge=False, max_depth=5000):
        """calculate cGENIE vertical depth
        :param N: number of grid points
        :param edge: if True, return edge points, otherwise return midpoints
        """
        ez0=0.1

        z1 = ez0 * ((1.0 + 1/ez0) ** (1.0 / N) - 1.0)
        dz = np.zeros(N)
        tv2 = 0.0

        for k in range(1, N + 1):
            tv3 = ez0 * ((z1 / ez0 + 1) ** k - 1)
            dz[N - k] = tv3 - tv2
            tv2 = tv3

        edges = np.zeros(N + 1)
        for k in range(N - 1, -1, -1):
            edges[k] = edges[k + 1] - dz[k]

        edges *= -max_depth  # Convert to physical depths (0 at surface, increasing positive down)

        if edge:
            return edges
        else:
            z = np.array([(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)])
            return z


    def get_normal_lon(self, N=36, edge=False):
        """
        Normal longitude in 10 degree resolution (default),
        if edge is False, then return midpoint
        """
        resolution = 360 / N
        if edge:
            lon_edge = np.linspace(-180, 180, N + 1)
            return lon_edge
        else:
            lon = np.linspace(-180+resolution/2, 180-resolution/2, N)
            return lon        
        
    def lon_n2g(self, x, grid_lon_offset=-260):
        """
        Convert normal longitude (-180, 180) to GENIE longitude (-270, 90)

        :param x: normal longitude
        :return: GENIE longitude
        """
        if grid_lon_offset == -180:
            print("this is already normal lat")
            return x
        normal_lon_cut = self.lon_g2n(grid_lon_offset)
        if x > normal_lon_cut and x < 180:
            return x - 360
        else:
            return x

    def lon_g2n(self, x):
        """
        Convert GENIE longitude to normal longitude.
        This is independent on the grid_offset_start option

        :param x: GENIE longitude
        :return: normal longitude        
        """

        if x < -180:
            return x + 360
        else:
            return x

    def lon_e2n(self, x):
        """ Convert eastern longitude to normal longitude

        :param x: longitude in eastern degree
        :return: normal longitude
        """
        
        if x > 180:
            return x - 360
        else:
            return x

    def lon_n2e(self, x):
        """Convert normal longitude (-180, 180) to longitude east(0,360)

        :param x: normal longitude
        :return: longitude in eastern degree
        """
        if x < 0:
            return 360 + x
        else:
            return x

    def xr_n2g(self, data: xr.Dataset, longitude="lon", *args, **kwargs) -> xr.Dataset:
        """Apply longitude conversion method n2g for the input data (normal to GENIE)

        :param data: input data
        :param longitude: longitude coordinate name

        :returns: xr.Dataset        
        """
        return data.assign_coords(
            {longitude: list(map(self.lon_n2g, data[longitude].values, *args, **kwargs))}
        ).sortby("lon")

    def xr_g2n(self, data: xr.Dataset, longitude="lon", *args, **kwargs) -> xr.Dataset:
        """Apply longitude conversion method g2n for the input data (GENIE to normal)

        :param data: input data
        :param longitude: longitude coordinate name

        :returns: xr.Dataset        
        """
        return data.assign_coords(
            {longitude: list(map(self.lon_g2n, data[longitude].values, *args, **kwargs))}
        ).sortby(longitude)

    def xr_e2n(self, data: xr.Dataset, longitude="lon", *args, **kwargs) -> xr.Dataset:
        """Apply longitude conversion method e2n for the input data (eastern to normal)

        :param data: input data
        :param longitude: longitude coordinate name
        
        :returns: xr.Dataset        
        """
        return data.assign_coords(
            {longitude: list(map(self.lon_e2n, data[longitude].values), *args, **kwargs)}
        ).sortby(longitude)

    def xr_n2e(self, data: xr.Dataset, longitude="lon", *args, **kwargs) -> xr.Dataset:
        """
        Apply longitude conversion method n2e for the input data (normal to eastern)

        :param data: input data
        :param longitude: longitude coordinate name

        :returns: xr.Dataset
        """
        return data.assign_coords(
            {longitude: list(map(self.lon_n2e, data[longitude].values, *args, **kwargs))}
        ).sortby(longitude)

    def mask_Arctic_Med(self, array, policy="na"):
        """
        mask Arctic and Meditterean Sea in cGENIE modern continent configuration

        :param array: 36x36 GENIE array
        """
        if policy == "na":
            array[34:36, :] = np.nan
            array[27:30, 25:30] = np.nan
        elif policy == "zero":
            array[34:36, :] = 0
            array[27:30, 25:30] = 0

        return array

    def GENIE_grid_mask(
        self, base="worjh2", basin="ALL", subbasin="", invert=False
    ):
        """
        Get a modern GENIE 36x36 mask array from input data.
        The input array is flipped (left/right flip -> up/down flip) for easy recognition

        :continent: worjh2, worlg4, worbe2, GIteiiaa, GIteiiaa, p0055c
        :basin: Atlantic/Pacific/Indian/ALL/Tanzania
        :subbasin: N/S/ALL, ALL means Southern Ocean section included

        :returns: GENIE grid array where continent/ice cap is 0 and ocean is 1, default is 'worjh2'
        """

        filename = f"mask_{base}_{basin}{subbasin}.txt"
        file_path=str(files('data').joinpath(filename))

        grid_mask_raw = np.loadtxt(file_path, dtype=int)
        grid_mask = np.flip(np.fliplr(grid_mask_raw))

        if invert:
            grid_mask = ~grid_mask + 2

        return grid_mask
        
    def geniebin_lat(self, x, *args,**kwargs):
        """
        Categorize <latitude> into cGENIE grid bins
        """
        ## check the latitude input range
        if x >= -90 and x <= 90:
            lat_edge = self.get_genie_lat(edge=True, *args,**kwargs)
            lat = self.get_genie_lat(edge=False)

            for i in range(36):
                if x > lat_edge[i] and x <= lat_edge[i + 1]:
                    x = lat[i]
        else:
            raise ValueError("Latitude must be in [-90,90]")

        return x

    def geniebin_lon(self, x, *args,**kwargs):
        """
        Categorize <longitude> into cGENIE grid bins
        """
        ## check the longitude input range
        if x >= -180 and x <= 180:
            lon_edge = self.get_genie_lon(edge=True, *args,**kwargs)
            if 'N' not in kwargs: N=36                
            for i in range(N):
                if x > lon_edge[i] and x <= lon_edge[i + 1]:
                    x = (lon_edge[i] + lon_edge[i + 1]) / 2  # middle value in the bin
        else:
            raise ValueError("Longitude must be in [-180,180]")

        return x

    def geniebin_depth(self, x, *args,**kwargs):
        """
        Categorize <depth> into cGENIE grid bins
        """
        ## check the depth input range
        if x >= 0 and x <= 5000:
            depth_edge = self.get_genie_depth(edge=True)
            depth = self.get_genie_depth(edge=False)

            for i in range(16):
                if x >= depth_edge[i] and x < depth_edge[i + 1]:
                    x = depth[i]
        else:
            raise ValueError("Depth must be in [0,5000]")

        return x

    def normbin_lon(self, x, *args,**kwargs):
        """
        Categorize <longitude> into normal grid bins
        """
        ## check the longitude input range
        if x >= -180 and x <= 180:
            lon_edge = self.get_normal_lon(edge=True, *args,**kwargs)
            if 'N' not in kwargs: N=36                
            for i in range(N):
                if x > lon_edge[i] and x <= lon_edge[i + 1]:
                    x = (lon_edge[i] + lon_edge[i + 1]) / 2
        else:
            raise ValueError("Longitude must be in [-180,180]")
        
        return x

    def haversine_distance(self, lat1, lon1, lat2, lon2):
        """
        Calculate the Haversine distance between corresponding pairs of points.

        :param lat1: Array of latitudes for points 1
        :param lon1: Array of longitudes for points 1
        :param lat2: Array of latitudes for points 2
        :param lon2: Array of longitudes for points 2
        
        :returns: Array of distances between corresponding pairs of points
        """
        # Convert latitude and longitude from degrees to radians
        lat1_rad, lon1_rad, lat2_rad, lon2_rad = np.radians(lat1), np.radians(lon1), np.radians(lat2), np.radians(lon2)

        # Haversine formula
        dlon = lon2_rad - lon1_rad
        dlat = lat2_rad - lat1_rad
        a = np.sin(dlat / 2) ** 2 + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2) ** 2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

        # Radius of the Earth in kilometers
        earth_radius = 6371.0
        distance = earth_radius * c

        return distance

    def geo_dis3d(self, point1, points2):
        """
        Calculate the 3D geographical distance between point1 and multiple points2.

        :param point1: Tuple/list of coordinates (z, lat, lon) or (lat, lon)
        :param points2: Numpy Array of shape (n, 3) containing coordinates (z, lat, lon) of points
        :returns: Array of distances between point1 and each point in points2
        """
        z1, lat1, lon1 = point1
        z2, lat2, lon2 = points2[:, 0], points2[:, 1], points2[:, 2]
        hor_dis = self.haversine_distance(lat1, lon1, lat2, lon2)
        ver_dis = np.abs(z1 - z2) / 1000  # Convert meters to kilometers
        return np.sqrt(hor_dis**2 + ver_dis**2)

    def geo_dis2d(self, point1, points2):
        """
        Calculate the 2D geographical distance between point1 and multiple points2.

        :param point1: Tuple/list of coordinates (z, lon, lat) or (lon, lat)
        :param points2: Array of shape (n, 3) containing coordinates (z, lon, lat) of points
        :returns: Array of distances between point1 and each point in points2
        """
        lat1, lon1 = point1[0], point1[1]
        lat2, lon2 = points2[:, 0], points2[:, 1]
        return self.haversine_distance(lat1, lon1, lat2, lon2)
    
    @staticmethod
    def check_dimension(input):
        """check the presence of latitude, longitude, depth, and time in the input tuple

        :param input: tuple of dimension names
        :return: tuple of boolean values indicating the presence of latitude, longitude, depth, and time

        Example
        -----------
        >>> input = ('lat', 'lon', 'depth', 'time')
        >>> check_dimension(input)
        (True, True, True, True)
        """
        # Convert tuple elements to lowercase
        input_lower = tuple(element.lower() for element in input)

        # Initialize flags for presence of each element
        has_lat = False
        has_lon = False
        has_depth = False
        has_time = False

        lat_candidates = ['lat', 'latitude', 'y', 'lat_moc', 'lat_psi', 'lat_opsi']
        lon_candidates = ['lon', 'longitude', 'x', 'lon_psi', 'lon_opsi']
        depth_candidates = ['depth', 'z', 'z_t', 'level', 'nlevel', 'lvl','lev', 'depth_1',
                            'elevation [m]', 'zt','depth_surface', 'zt_moc', 'deptht']
        time_candidates = ['time', 't', 'age', 'date', 'year', 'age [ka]', 'age [ka bp]']

        # Check for presence of each element
        for element in input_lower:
            if element in lat_candidates:
                has_lat = True
            if element in lon_candidates:
                has_lon = True
            if element in depth_candidates:
                has_depth = True
            if element in time_candidates:
                has_time = True

        return has_lat, has_lon, has_depth, has_time

    def dim_order(self,input):
        """
        Determine the order of dimensions in the input data array

        :return: tuple of index in the input data

        Example
        -----------
        >>> input = ('lat', 'lon', 'depth')
        >>> dim_order(input) ## always follow time, depth, lat, lon
        (2, 0, 1) 
        """
        order = []
        input_lower = tuple(element.lower() for element in input)

        dim_candidates = [
            ['time', 't', 'age', 'date', 'year', 'age [ka]', 'age [ka bp]'],
            ['depth', 'z', 'z_t', 'level', 'nlevel', 'lvl', 'lev', 'depth_1',
                'elevation [m]','zt', 'depth_surface','zt_moc'],
            ['lat', 'latitude', 'y', 'lat_moc', 'lat_psi'],
            ['lon', 'longitude', 'x', 'lon_psi']
        ]

        for candidates in dim_candidates:
            for candidate in candidates:
                if candidate in input_lower:
                    order.append(input_lower.index(candidate))
                    break

        return tuple(order)
    
    @staticmethod
    def set_coordinates(obj, index):
        """Set the coordinates of the input object based on the input index

        :param obj: object to set coordinates
        :param index: tuple of dimension names        
        """
        obj.n_index = len(index)
        gp = GridOperation()
        has_lat, has_lon, has_depth, has_time = gp.check_dimension(index)
        index_order = gp.dim_order(index)
        match obj.n_index:
            case 1:
                if has_lat: obj.lat = index[0]
                if has_lon: obj.lon = index[0]                
                if has_depth: obj.depth = index[0]
                if has_time: obj.time = index[0]
            case 2:
                if has_lat and has_lon:
                    obj.lat = index[index_order[0]]
                    obj.lon = index[index_order[1]]
                if has_lat and has_depth:
                    obj.depth = index[index_order[0]]
                    obj.lat = index[index_order[1]]
                if has_lat and has_time:
                    obj.time = index[index_order[0]]
                    obj.lat = index[index_order[1]]
                if has_lon and has_depth:
                    obj.depth = index[index_order[0]]
                    obj.lon = index[index_order[1]]
                if has_lon and has_time:
                    obj.time = index[index_order[0]]
                    obj.lon = index[index_order[1]]
                if has_depth and has_time:
                    obj.time = index[index_order[0]]
                    obj.depth = index[index_order[1]]
            case 3:
                if not has_lat:
                    obj.time = index[index_order[0]]
                    obj.depth = index[index_order[1]]
                    obj.time = index[index_order[0]]                                          
                if not has_lon:
                    obj.time = index[index_order[0]]
                    obj.depth = index[index_order[1]]
                    obj.lat = index[index_order[2]]
                if not has_depth:
                    obj.time = index[index_order[0]]
                    obj.lat = index[index_order[1]]
                    obj.lon = index[index_order[2]]
                if not has_time:
                    obj.lat = index[index_order[0]]
                    obj.lon = index[index_order[1]]
                    obj.depth = index[index_order[2]]
            case 4:
                    obj.time = index[index_order[0]]
                    obj.depth = index[index_order[1]]
                    obj.lat = index[index_order[2]]
                    obj.lon = index[index_order[3]]             

class Interpolator:

    """A univeral ineteprolator for cgeniepy that
    can be used to interpolate data for both regular and irregular grid
    """
    
    def __init__(self, dims, coordinates, values, grid_number=200, method="r-linear"):
        """
        initialize regridder

        :param dims: dimensions of the data
        :param coordinates: a list of coordinate arrays
        :param values: data values
        :param grid_number: number of grid points
        :param method: interpolation method, use "x-y" format. where x is strcutre or not, y is the algorithm
        """
        self.dims = dims
        self.coords = coordinates
        self.values = values

        ## number of grid points
        self.grid_number = grid_number
        ## create new coordinates
        self.gridded_coord = self.new_coordinate(n=grid_number)
        ## create meshgrid
        self.meshgrid = self.new_meshgrid(*self.gridded_coord, indexing="ij")
        ## interpolation function
        self.interp_function = self._create_interp_function(method=method)
        ## interpolate data
        self.gridded_data = self.interpolate_data(tuple(self.meshgrid))

    def _create_interp_function(self, method):
        """
        create interpolation function

        :param method: interpolation method, use "x-y" format. where x is strcutre or not, y is the algorithm
        """

        if method=='1d':
            ## use cubic spline
            return CubicSpline(self.coords[0], self.values)

        data_class = method.split("-")[0]
        ## regular grid interpolation
        if data_class == "r":
            true_method = method.split("-")[1]
            interp_function = RegularGridInterpolator(
                self.coords, self.values, method=true_method
            )
        ## irregular grid
        elif data_class == "ir":
            true_method = method.split("-")[1]
            if true_method == "linear":
                interp_function = LinearNDInterpolator(self.coords, self.values)
            elif true_method == "nearest":
                interp_function = NearestNDInterpolator(self.coords, self.values)
        else:
            raise ValueError("Method not supported")

        return interp_function

    def new_coordinate(self, n):
        """create new coordinates for regridding
        
        :param n: the resolution of the new grid

        :returns: a list of new coordinates
        """
        new_coords = []
        for coord_values in self.coords:
            min_val, max_val = np.min(coord_values), np.max(coord_values)
            new_values = np.linspace(min_val, max_val, n)
            new_coords.append(new_values)
        return new_coords

    def new_meshgrid(self, *args, **kwargs):
        """create meshgrid for regridding

        :returns: numpy meshgrid
        """

        return np.meshgrid(*args, **kwargs)

    def interpolate_data(self, *args, **kwargs):
        """Use the interpolation function to get regridded values

        :returns: regridded data array
        """
        return self.interp_function(*args, **kwargs)

    def to_xarray(self):
        """
        Combine the regridded data with the new coordinates to create a xarray DataArray

        :returns: xarray DataArray
        """
        return xr.DataArray(
            self.gridded_data, dims=self.dims, coords=self.gridded_coord
        )

    def to_dataframe(self):
        """
        Convert the regridded data to a pandas DataFrame

        :returns: pandas DataFrame
        """
        assert len(self.dims) == 1, "Only support 1D data"

        index = self.gridded_coord[0]
        values = self.gridded_data[0]

        ## build a data frame
        data = {self.dims[0]: index, "interpolated_values": values}

        return pd.DataFrame(data)
    
