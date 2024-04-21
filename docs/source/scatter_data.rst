Scatter Data
===================


Basic concept
-------------------
Basides the GriddedData that stores array-like data, the ScatterData class is used to store scattered data in the form of data frame. It is a container of pandas DataFrame object that stores the data in a table format, with some optimsation for coordinate-based operations. The ScatterData and GriddedData classes can be transformed to each other.



Read the data into ScatterData
--------------------------------

ScatterData can be created from a pandas DataFrame object or directly from a file. It accepts normal csv and txt file, and the tab file as one downloads from Pangeae. 
The only difference from pandas DataFrame is that the ScatterData functionalities are based on their index columns. So one must set the index columns before doing further operations.

.. code-block:: python

    from cgeniepy.table import ScatterData

    sd = ScatterData('data.csv')
    sd.set_index(['lat', 'lon', 'depth'])


Detect basin based on coordinates
-----------------------------------
A fun method to use is the detect_basin method, which can detect the basin based on the coordinates. It is useful when one wants to know the basin of a certain point. The method returns the basin name and the basin code.

.. code-block:: python

    sd.detect_basin() ## it will automatically detect the basin based on the coordinates and add the basin column into the data frame