Extract simulated time-slice and time-series data
========================================================

Once initialised the GenieModel object, it will use the model path and gemflag to find netcdf/timeseries file. Then it stores all the available variables and corresponding path in the attribute `nc_vardict` and `ts_varlist`. By such, users can forget about the model path and focus on the target variable to analyse.

The model instance provides two methods to access the model time-slice and time-series respectively: `get_var`and `get_ts`. Each receives a parameter of the variable name. The variable name here means the short name in the netcdf file. For example, `ocn_temp` stands for ocean temperature.

Read time slice data
-------------------------------
.. code-block:: python

    model.get_var("ocn_temp")

		
Read time series data
-------------------------------
.. code-block:: python

    model.get_ts("ocn_temp")

Read diagnostic table
-------------------------------
It is worth noting that cGENIE also produces a global average diagnostic table. This table contains the global average of all the useful variables in the model. One can use `get_diag_avg` to read this information.

.. code-block:: python

    model.get_diag_avg(9999) ## get the summary of the 9999th year
