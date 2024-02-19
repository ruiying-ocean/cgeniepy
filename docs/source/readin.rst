Read in model output
====================

Download from cluster
---------------------

The first step starts from getting data from the cluster. If it is in `tar.gz` format then uncompress it. Otherwise, directly download as a folder.

.. code_block::

   $ scp -r remote:cgenie_output/experiment_id ~/Downloads
   

Initialise a model instance
-------------------------------

.. code-block:: python

    from cgeniepy.model import GenieModel

    # Single model
    model = GenieModel('path_a')

    # Model ensemble
    multi_dirs = ['/path_a/', '/path_b/']
    model = GenieModel(path_to_model_output)

Read time slice data
-------------------------------

.. code-block:: python

    model.get_var("ocn_temp").array

		
Read time series data
-------------------------------
.. code-block:: python

    model.get_var("biogem_series_ocn_temp.res").array
