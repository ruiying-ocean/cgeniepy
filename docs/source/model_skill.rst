Model-data comparison
==============================


ArrComparison class
-----------------------
The ArrComparison class is provided to compare the model and data when both are in the form of arrays.

.. code-block:: python

    from cgeniepy.skill imoprt ArrComparison, DFComparison
    ac = ArrComparison(model_array, data_array)
    ac.mscore() ## M-score
    ac.rmse()  ## RMSE
    ac.pearson_r() ## Pearson's correlation coefficient
    ac.plot() ## plot the model and data arrays in 1:1 line plot



DFComparison class
-----------------------
The DFComparison class is provided to compare the model and data when they are the different columns in the data frame. Other methods are same to ArrComparison class because they are parent-child classes.

.. code-block:: python

    dfc = DFComparison(dataframe, model_column, data_column)
    dfc.mscore() ## M-score
    dfc.rmse()  ## RMSE
    dfc.pearson_r() ## Pearson's correlation coefficient
    dfc.plot() ## plot the model and data arrays in 1:1 line plot
