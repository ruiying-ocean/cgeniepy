Installation
==============

.. _installation:

Install Python
----------------
This is a Python package, so you need to have Python in your computer. As the most popular programming language, I assume most people can find out how to install Python by simply googling it. Visiting the `Python official website <https://www.python.org/downloads/>`_ should be enough to download the latest Python (> 3.10). Make sure it is included in the shell environment (which means you can type :code:`python` in your terminal without triggering a error).

Use pip
----------------

`pip` is an official package-management system in Python ecosystem. Once you have download Python, you can use `pip` to download millions of open-source packages. Usually, you just need to type three words to install one package:

Open *terminal* (MacOS, GNU/Linux) or `cmd` (Win) and type:

.. code-block:: console

       $ pip install package_name

Install cgeniepy
--------------------

Now let's install `cgeniepy` using `pip`. It is slightly different from the last command because I want to avoid any possible error in the installation step. This command will download the package that I built and uploaded to testpypi using the specific Python3 version (in some old operation system python points to python2).

.. code-block:: console

       $ python3 -m pip install -i https://test.pypi.org/simple/ cgeniepy


Alternatively, we can directly download the developing version from GitHub:

.. code-block:: console
		
       $ python3 -m pip install git+https://github.com/ruiying-ocean/cgeniepy.git@master


(Optional) Install dependencies       
---------------------------------
cgeniepy is not a standalone software. It is built on other excellent packages (so called "dependencies"). Although installing cgeniepy will automatically download the related dependent packages, sometime it could fail. Using the blow command can manually avoid this.

.. code-block:: console

       $ python3 -m pip install xarray numpy matplotlib pandas scipy

This also touches the most famous Python data science-related packages:

- **xarray**: N-D labeled arrays and datasets in Python, I use it to handle the NetCDF data.
- **matplotlib**: A comprehensive library for creating static, animated, and interactive visualizations in Python. This is the plotting engine in cgeniepy.
- **numpy**: The fundamental package for scientific computing with Python. 
- **pandas**: A fast, powerful, flexible and easy to use open source data analysis and data manipulation library built on top of Python. It is the core of cgeniepy.ScatterData and used to deal with the cGENIE timeseries output.
- **scipy**: A Python-based ecosystem of open-source software for mathematics, science, and engineering.


(Optional) Verify the installation
-------------------------------------
One can use pytest to do the unit test. If you have pytest installed, you can run the following command to test the installation:

.. code-block:: console

       $ pytest --pyargs cgeniepy
