Computational analysis of time-slice data
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


Selecting data
-----------------------
The selection of data can be done by using the `sel` method or `isel` method. The `sel` method is used to select the data by the coordinate value, while the `isel` method is used to select the data by the index. Note that the selection and any computation of data can be inplace or not. If you want to keep the original data, you need to change the attribute `mutable` to `False`.


.. code-block:: python

    sst = model.get_var("ocn_sur_temp", mutable=False)
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
Similar to the selection of data by coordinate (time, lat, long etc), you can mask a ocean basin in cgeniepy. The first method is to use `mask_basin` method. It reads the pre-stored basin mask for the specific basic configration (e.g., 'worjh2' and 'worlg4' for modern ocean topography in cGNEIE).


.. code-block:: python

    sst = model.get_var("ocn_sur_temp", mutable=False)
    sst.mask_basin(base="worjh2", basin='Atlantic') ## -> mask the other oceans except Atlantic basin


The other way is to use `sel_modern_basin` method. As the name suggests, it only works for the modern model output. In fact, it is based on the basin divsion in IPCC AR6 and the provided functionalisties in `regionmask` package. The only caveat is that it only works for lat-lon data.


.. code-block:: python

    sst = model.get_var("ocn_sur_temp", mutable=False)
    sst.sel_modern_basin('NPO') ## -> select the North Pacific Ocean


Chain computation
-----------------------
All the methods can be done in a chain. For example, you can select the data, calculate the mean value and plot it in a single line. The only thing to remember is that it change the data in place if the attribute `mutable` is `True`.


.. code-block:: python

    sst = model.get_var("ocn_sur_temp")
    sst.sel_modern_basin('NPO').mean() ## -> select the data in the northern hemisphere, calculate the mean value

