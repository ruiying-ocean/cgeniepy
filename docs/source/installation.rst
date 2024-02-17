Installation
============

.. _installation:

Install Python
--------------
First visit `Python website <https://www.python.org/downloads/>`_ to download the latest Python (> 3.10). Make sure it is included in the environment (echo $PATH).

Install cgeniepy
----------------

To use cgeniepy, first install it using pip:

.. code-block:: console

       $ python3 -m pip install -i https://test.pypi.org/simple/ cgeniepy


Alternatively, we can directly download the developing version from GitHub:

.. code-block:: console
		
       $ python3 -m pip install git+https://github.com/ruiying-ocean/cgeniepy.git@master


Install dependencies       
--------------------
.. code-block:: console

       $ python3 -m pip install xarray numpy matplotlib pandas scipy
