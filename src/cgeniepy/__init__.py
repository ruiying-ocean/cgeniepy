# from pint import UnitRegistry
# import pathlib


# ureg = UnitRegistry()
# Q_ = ureg.Quantity
# file_path = pathlib.Path(__file__).parent.parent / "data/context.txt"
# ureg.load_definitions(file_path)


from importlib.resources import files
from cgeniepy.model import GenieModel
from cgeniepy.ecology import EcoModel
from cgeniepy.plot import CommunityPalette

def sample_model(model_type='GenieModel', *args, **kwargs):
    file_path=str(files('data').joinpath('sample_model'))
    
    if model_type == 'GenieModel':
        model = GenieModel(file_path,*args, **kwargs)
    elif model_type == 'EcoModel':
        model = EcoModel(file_path,*args, **kwargs)
    else:
        raise ValueError('model_type must be either GenieModel or EcoModel')    

    return model

def load_obs(grid='worjh2'):
    file_path=str(files('data').joinpath(grid+'_obs.nc'))
    import xarray as xr
    obs_data= xr.load_dataset(file_path)
    return obs_data

def register_cmap():
    import matplotlib as mpl

    available_cmaps = CommunityPalette().avail_palettes()
    for cmap_name in available_cmaps:
        community_cmap = CommunityPalette(name=cmap_name).colormap
        mpl.colormaps.register(name=cmap_name, cmap=community_cmap)

try:
    register_cmap()
except Exception as e:
    print(f"Error registering colormaps: {e}")
