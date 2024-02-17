Usage
=====

Download model output
---------------------

The first step starts from getting data from the cluster. If it is in `tar.gz` format then uncompress it. Otherwise, directly download as a folder.

.. code_block::

   $ scp -r remote:cgenie_output/experiment_id ~/Downloads
   

Explore the model output
------------------------

Firstly import the package

.. code_block:: Python

   $ from cgeniepy.model import ModelOutput
   
