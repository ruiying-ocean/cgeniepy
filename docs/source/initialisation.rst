Initialise a model experiment instance
===========================================

Download cGENIE output from cluster
--------------------------------------

The first step starts from getting data from the cluster. If it is in `tar.gz` format then uncompress it. Otherwise, directly download as a folder.

.. code-block:: console

   $ scp -r remote:cgenie_output/experiment_id ~/Downloads

For test, an example output is already attached with the package under directory `cgeniepy/src/data/muffin.CBE.worlg4.BASESFeTDTL.SPIN`. This is an modern configuration with ecosystem component enabled.


Initialise a model instance
-------------------------------
The next step is to initialise a model experiment instance. Two parameters are needed: (1) the directory path (e.g., `muffin.CBE.worlg4.BASESFeTDTL.SPIN`); (2) the target sub-model (e.g., biogem or ecogem). The default of sub-model is biogem, but if you want to access more data, you can create of gemflag.

.. code-block:: python

    from cgeniepy.model import GenieModel
    
    model = GenieModel('/Users/XX/Downloads/path_a', gemflag=['biogem', 'ecogem'])

The model also support multiple experiment directory. The only change is to use a list of directories as input.
This list can be create by yourself, or use Python's Pathlib module.

.. code-block:: python

    # Model ensemble    
    model_dir = ['/path_a/', '/path_b/']
    
    models = GenieModel(model_dir)

    ## Alternatively
    from pathlib import Path
    model_dir = Path("XX/model")
    models = GenieModel(model_path)


Read a remote OPeNDAP dataset
-----------------------------

Pass the HTTP or HTTPS URL of an OPeNDAP dataset directly to ``GenieModel``.
Use the data endpoint, rather than a catalogue page, an HTML preview, or an
experiment directory. The endpoint does not need to end in ``.nc``.

.. code-block:: python

    from cgeniepy import GenieModel

    url = "https://your-server/thredds/dodsC/run/fields_biogem_2d.nc"
    with GenieModel(url) as model:
        print(model.ncvar_dict)  # Discover the available variables.
        temperature = model.get_var("ocn_sur_temp").isel(time=-1)
        temperature.data.load()

Opening a single dataset reads its metadata and coordinates. Field values are
read as needed, so select the time, depth, or region before loading or computing.
Call ``.data.load()`` inside the ``with`` block if you want to keep using the
selected data after the connection closes.

To try this with a public dataset, NOAA provides
`ERSST v5 sea-surface temperatures <https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html>`_:

.. code-block:: python

    url = "https://psl.noaa.gov/thredds/dodsC/Datasets/noaa.ersst.v5/sst.mnmean.nc"
    with GenieModel(url) as model:
        sst = model.get_var("sst").isel(time=-1)
        pacific = sst.sel(lat=0, lon=180, method="nearest")
        pacific.data.load()
        print(pacific.data)  # Temperature in degrees Celsius at this grid point.

The default remote reader is ``netcdf4``, which is already a cgeniepy dependency.
Options such as ``decode_times=False`` are passed to ``xarray.open_dataset``.
For a different reader or an authenticated session, install ``pydap`` and pass
the corresponding xarray options:

.. code-block:: python

    # Install the optional reader with: pip install pydap
    with GenieModel(url, engine="pydap", decode_times=False) as model:
        temperature = model.get_var("sst").isel(time=-1)
        temperature.data.load()

For authentication, an existing authenticated ``requests.Session`` can be passed
as ``backend_kwargs={"session": session}`` with ``engine="pydap"``. See the
`xarray OPeNDAP guide <https://docs.xarray.dev/en/latest/user-guide/io.html#opendap>`_
for reader and authentication details.

A list or tuple of dataset URLs creates an ensemble along the ``model``
dimension, just like a list of local directories. Each URL represents a separate
model with matching variables; a list is not a way to combine the 2-D and 3-D
files of one experiment. Ensemble concatenation can load field values; install
``dask[array]`` and pass ``chunks={}`` to keep ensemble arrays lazy. Local
directories and remote URLs cannot be mixed in one ensemble.

Remote access covers the variables in the supplied dataset. ``gemflag`` only
controls discovery in local directories. ``get_ts()``, ``get_diag_avg()``, and
``get_config()`` require local text files. ``EcoModel`` also requires a local
``Plankton_params.txt``; use ``GenieModel.get_var()`` to read remote ECOGEM fields.
