from matplotlib.testing.decorators import image_comparison
from cgeniepy.table import ScatterData
from importlib.resources import files
import cgeniepy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def create_sample_data():
    model = cgeniepy.sample_model()
    return model.get_var('ocn_sur_temp').isel(time=-1)

@image_comparison(baseline_images=['test_map'], remove_text=True,
                  extensions=['png'], style='mpl20')
def test_map():    
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Mollweide()})
    data = create_sample_data()
    data.plot(ax=ax)
    
    return fig

@image_comparison(baseline_images=['test_line'], remove_text=True,
                  extensions=['png'], style='mpl20')
def test_line():
    data = create_sample_data()
    fig, ax = plt.subplots()
    data.mean(dim='lon').plot(ax=ax)
    return fig

@image_comparison(baseline_images=['test_scatterdatavis'], remove_text=True,
                  extensions=['png'], style='mpl20')
def test_scatterdatavis():
    file_path = str(files("data").joinpath("EDC_CO2.tab"))
    data= ScatterData(file_path, sep='\t')
    data.set_index('Age [ka BP]')
    fig, ax = plt.subplots()
    data.plot(var='CO2 [ppmv]', ax=ax)
    return fig