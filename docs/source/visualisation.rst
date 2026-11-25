Visualisation of time-slice data
=====================================


Basic concept
-----------------
The GriddedData has a method `plot` to visualise the data. It will calculate the data dimension and plot the data accordingly. For example, if it is a 2D data with latitude and longitude, it will plot the data on a map. If it is a 3D data with latitude, longitude and time, it will plot many maps for each time slice. The plotting engine is `matplotlib` and `cartopy`.


Module-based plotting design
-------------------------------
The plot method is designed to be modular. You can choose different element: contour, pcolormesh (default), colorbar, land outline, gridline etc.

With ``contourf=True``, the filled contours reach the model's coastline and the land cells are drawn over them, so pcolormesh is off unless you pass ``pcolormesh=True``. Contour lines drawn with them follow the filled levels, unless you set ``contour_kwargs["levels"]``.


.. code-block:: python

    sst = model.get_var("ocn_temp")
    sst.plot(colorbar=True, outline=True)


Customisation
----------------
I have to admit it is not very smart to do the customisation. Instead of like `matplotlib` which you can pass all the arguments to the plotting function, I store the parameters in a dictionary to help achieve the modular design. For example, if you want to change the color map of pcolormesh plot, do it like this:

.. code-block:: python

    p = sst.to_GriddedDataVis()
    p["pcolormesh_kwargs"]["cmap"] = "Spectral_r"
    p.plot()



Map frame in the style of NCL and GMT
---------------------------------------
On a rectangular projection (``PlateCarree``, ``Miller``, ``Mercator``), ``gridline=True`` draws the frame that NCL and GMT draw: tick marks pointing out on all four sides, labelled in degrees on the bottom and the left, with two minor ticks between the major ones, over dotted grid lines. ``zebra_frame=True`` adds GMT's fancy black-and-white frame on top of those ticks. Curved projections such as the default Equal Earth get the grid lines only.

.. code-block:: python

    import cartopy.crs as ccrs

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})
    sst.plot(ax=ax, contourf=True, colorbar=True, gridline=True)        # NCL
    sst.plot(ax=ax, colorbar=True, gridline=True, zebra_frame=True)     # GMT

The tick spacing comes from ``gridline_kwargs["xlocs"]`` and ``gridline_kwargs["ylocs"]`` (every 60° and 30° by default), and ``gridline_kwargs["draw_labels"] = False`` keeps the ticks without their labels. ``cgeniepy.plot.add_border_ticks(ax)`` draws the same frame on any rectangular cartopy map.

GENIE's own grid view
---------------------
``projection="sine"`` draws the map as GENIE's equal-area grid is built, and as genieplot draws it: latitude spaced by its sine, so every row of the 36 x 36 grid is equally tall and every cell takes the same area on the page. The map is centred on the grid's own longitudes (-260 to 100 °E, so centred on 80 °W), in a 2:1 box. It is cartopy's cylindrical equal-area projection (``LambertCylindrical``) on a sphere, so the frame labels of ``gridline=True`` and ``zebra_frame=True`` work on it. ``projection`` also takes any cartopy CRS for the map's new axes; with ``ax`` the axes bring their own.

.. code-block:: python

    sst.plot(projection="sine", colorbar=True, gridline=True)

Near the poles the latitude labels crowd together (60° is at 0.87 of the height); ``p.aes_dict["gridline_kwargs"]["ylocs"] = [-60, -30, 0, 30, 60]`` leaves out the poles.

More Examples
------------------
Please see the gallery section.
