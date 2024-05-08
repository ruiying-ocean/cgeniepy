Initialise a model experiment instance
===========================================

Download cGENIE output from cluster
--------------------------------------

The first step starts from getting data from the cluster. If it is in `tar.gz` format then uncompress it. Otherwise, directly download as a folder.

.. code-block:: console

   $ scp -r remote:cgenie_output/experiment_id ~/Downloads
   

Initialise a model instance
-------------------------------
The next step is to initialise a model experiment instance. Two parameters are needed: (1) the directory path; (2) the target sub-model (e.g., biogem or ecogem).
The default of sub-model is biogem, but if you want to access more data, you can create of gemflag.

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
