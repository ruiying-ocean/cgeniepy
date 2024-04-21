visualisation of time-slice data
=====================================


Basic concept
-----------------
The GriddedData has a method `plot` to visualise the data. It will calculate the data dimension and plot the data accordingly. For example, if it is a 2D data with latitude and longitude, it will plot the data on a map. If it is a 3D data with latitude, longitude and time, it will plot many maps for each time slice. The plotting engine is `matplotlib` and `cartopy`.


Module-based plotting design
-------------------------------
The plot method is designed to be modular. You can choose different element: contour, pcolormesh (default), colorbar, land outline, gridline etc.


.. code-block:: python

    sst = model.get_var("ocn_temp")
    sst.plot(colorbar=True, outline=True)



More Examples
-----------
Please see the gallery section.