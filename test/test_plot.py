import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cgeniepy.table import ScatterData
import cgeniepy

def create_sample_data():
    model = cgeniepy.sample_model()
    return model.get_var('ocn_sur_temp').isel(time=-1)

def test_map_creation():
    """Test that map plotting doesn't crash and creates expected elements"""
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Mollweide()})
    data = create_sample_data()
    
    # Test that plotting works without error
    data.plot(ax=ax)
    
    # Test that the plot has expected properties
    assert len(ax.collections) > 0  # Has plot elements
    
    plt.close(fig)

def test_line_plot():
    """Test line plotting functionality"""
    data = create_sample_data()
    fig, ax = plt.subplots()
    
    line_data = data.mean(dim='lon')
    line_data.plot(ax=ax)
    
    # Test plot properties
    assert len(ax.lines) > 0  # Has line elements
    assert ax.get_xlabel() != ""  # Has x-label
    assert ax.get_ylabel() != ""  # Has y-label
    
    plt.close(fig)

def test_scatterdata_plot():
    """Test scatter data visualization"""
    
    from importlib.resources import files
    ## access the sample data file
    ## but convert to str for the ScatterData class
    file_path = str(files('cgeniepy.data').joinpath('EDC_CO2.tab'))

    data = ScatterData(file_path, sep='\t')
    data.set_index('Age [ka BP]')
    
    fig, ax = plt.subplots()
    data.plot(var='CO2 [ppmv]', ax=ax)
    
    # Test that data was plotted
    assert len(ax.lines) > 0 or len(ax.collections) > 0
    
    plt.close(fig)
