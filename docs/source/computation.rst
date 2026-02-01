Basic computation of time-slice data
===============================================================

Basic concept
--------------
The timeslice data in the model output is stored as GriddedData object. This is in essential a container of xarray DataArray object by storing it in the attribute `data`. If you don't know xarray, you can think of it as a N-dimensional array with additional metadata such as unit, long name. However, by design GriddedData provides additional functionalisties to manipulate the data, such as prettier plot, finding the nearest point, etc.


.. code-block:: python

    sst = model.get_var("ocn_sur_temp")
    sst.data ## -> a xarray DataArray object


Basic computation
-----------------------
Any basic computation can be done like a normal xarray DataArray object. For example, you can calculate the mean, maximum, minimum, standard deviation, etc.


.. code-block:: python

    sst = model.get_var("ocn_sur_temp")
    sst + 273.15 ## -> convert to Kelvin
    sst.mean() ## -> calculate the mean value
    sst.max() ## -> calculate the maximum value
    sst.min() ## -> calculate the minimum value
    sst.sd() ## -> calculate the standard deviation

Weighted mean
-----------------------

.. code-block:: python
		
    model = GenieModel("xxx", gemflag='biogem')

    o2 = model.get_var('ocn_O2').isel(time=-1) ##mol/kg
    ocn_vol = model.grid_volume().isel(time=-1) ##m3

    print("average of o2 weighted by ocean grid volume", o2.weighted_mean(ocn_vol.data.values))

    ## unweighted average of o2
    print("unweighted average of o2", o2.data.mean().values)

Selecting data
-----------------------
The selection of data can be done by using the `sel` method or `isel` method. The `sel` method is used to select the data by the coordinate value, while the `isel` method is used to select the data by the index. Note that the selection and any computation of data can be inplace or not.


.. code-block:: python

    sst = model.get_var("ocn_sur_temp")
    sst.isel(time=-1) ## -> select the last time slice
    sst.sel(sst.data.lat > 0) ## -> select the data in the northern hemisphere


Search the nearest point
----------------------------
It is useful to do the model-data comparison by finding the nearest point in the model output. The `search_point` method is designed to achieve this. By default, it uses xarray's `sel` method to find the nearest point. However, by passing the argument `ignore_na=True`, it will ignore the missing value in the data. This is my own implementation (inspired by the issue in xarray) based on geodistance calculation. But of course it is significantly slower than the xarray's method.

The only input is the coordinate of the point you want to search in the order of data's dimension. For example, if the data is 3D (time, lat, lon), you need to pass the coordinate in the order of (time, lat, lon).

.. code-block:: python
    
    sst = model.get_var("ocn_sur_temp")
    point = (0, 50) ## lat, lon
    sst.search_point(point)


Mask data
-----------------------
Similar to the selection of data by coordinate (time, lat, long etc), you can mask a ocean basin in cgeniepy.

Method 1: Using pre-defined cGENIE basin masks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The first method is to use `mask_basin` method. It reads the pre-stored basin mask for the specific basic configration (e.g., 'worjh2' and 'worlg4' for modern ocean topography in cGENIE).


.. code-block:: python

    sst = model.get_var("ocn_sur_temp")
    sst.mask_basin(base="worjh2", basin='Atlantic') ## -> mask the other oceans except Atlantic basin

    ## You can also combine multiple basins
    sst.mask_basin(base="worjh2", basin=['Atlantic', 'Pacific'])

Method 2: Using IPCC AR6 basin definitions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The other way is to use `sel_modern_basin` method. As the name suggests, it only works for the modern model output. In fact, it is based on the basin division in IPCC AR6 and the provided functionalities in `regionmask` package. The only caveat is that it only works for lat-lon data.

**Available basins:**

The AR6 ocean basins include:

- **46**: Arctic Ocean (AO)
- **47**: North Pacific Ocean (NPO)
- **48**: Equatorial Pacific Ocean (EPO)
- **49**: South Pacific Ocean (SPO)
- **50**: North Atlantic Ocean (NAO)
- **51**: Equatorial Atlantic Ocean (EAO)
- **52**: Southern Atlantic Ocean (SAO)
- **53**: North Indian Ocean (NIO)
- **55**: Equatorial Indian Ocean (EIO)
- **56**: South Indian Ocean (SIO)
- **57**: Southern Ocean (SO)

You can use either basin indices (int) or abbreviations (str):

.. code-block:: python

    sst = model.get_var("ocn_sur_temp")

    ## Using basin index
    sst.sel_modern_basin(47) ## -> select the North Pacific Ocean

    ## Using basin abbreviation
    sst.sel_modern_basin('NPO') ## -> select the North Pacific Ocean

    ## Select multiple basins to combine regions
    sst.sel_modern_basin(['NAO', 'EAO', 'SAO']) ## -> entire Atlantic Ocean
    sst.sel_modern_basin([47, 48, 49]) ## -> entire Pacific Ocean

    ## You can also chain with other operations
    atlantic_mean_sst = sst.sel_modern_basin([50, 51, 52]).mean()

**Reference:** The basin definitions follow Iturbide et al., (2020) ESSD and can be visualized at: https://regionmask.readthedocs.io/en/stable/_images/plotting_ar6_all.png

**See also:** :ref:`plot_basin_detection` example for comprehensive visualization of basin detection features.


Chain computation
-----------------------
All the methods can be done in a chain. For example, you can select the data, calculate the mean value and plot it in a single line. 


.. code-block:: python

    sst = model.get_var("ocn_sur_temp")
    sst.sel_modern_basin('NPO').mean() ## -> select the data in the northern hemisphere, calculate the mean value

