
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "auto_examples/plot_logo.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        :ref:`Go to the end <sphx_glr_download_auto_examples_plot_logo.py>`
        to download the full example code.

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_plot_logo.py:


================================
Customise the 2D map projection
================================

This example shows how to customise the 2D map including the projection, the color map, which is used as the logo of this package.

.. GENERATED FROM PYTHON SOURCE LINES 8-26



.. image-sg:: /auto_examples/images/sphx_glr_plot_logo_001.png
   :alt: plot logo
   :srcset: /auto_examples/images/sphx_glr_plot_logo_001.png
   :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 .. code-block:: none

    /Users/yingrui/cgeniepy/src/cgeniepy/model.py:51: UserWarning: No gemflag is provided, use default gemflags: [biogem]
      warnings.warn("No gemflag is provided, use default gemflags: [biogem]")






|

.. code-block:: Python

    import cgeniepy
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    ## Read in the model
    model = cgeniepy.sample_model()
    sst = model.get_var("ocn_sur_temp").isel(time=-1)

    ## use the Orthographic projection
    ## a full list of projections can be found at
    ## https://scitools.org.uk/cartopy/docs/latest/reference/projections.html#cartopy-projections
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Orthographic()})

    ## set the color map
    sst_plotter= sst.to_GriddedDataVis()
    sst_plotter.aes_dict['pcolormesh_kwargs']['cmap'] = plt.cm.inferno
    sst_plotter.plot(ax=ax, outline=True)
    plt.show()


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** (0 minutes 0.188 seconds)


.. _sphx_glr_download_auto_examples_plot_logo.py:

.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-example

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download Jupyter notebook: plot_logo.ipynb <plot_logo.ipynb>`

    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download Python source code: plot_logo.py <plot_logo.py>`

    .. container:: sphx-glr-download sphx-glr-download-zip

      :download:`Download zipped: plot_logo.zip <plot_logo.zip>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
