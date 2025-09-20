from importlib.resources import files
from . import data as _data_pkg
from cgeniepy.model import GenieModel
from cgeniepy.ecology import EcoModel
from cgeniepy.plot import CommunityPalette

def sample_model(model_type='GenieModel', *args, **kwargs):
    # Resolve from the actual data package module for robustness
    file_path = files(_data_pkg).joinpath('sample_model')
    
    
    if model_type == 'GenieModel':
        model = GenieModel(file_path,*args, **kwargs)
    elif model_type == 'EcoModel':
        model = EcoModel(file_path,*args, **kwargs)
    else:
        raise ValueError('model_type must be either GenieModel or EcoModel')    

    return model

def load_obs(grid='worjh2'):
    file_path = files(_data_pkg).joinpath(grid+'_obs.nc')
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
except ImportError as e:
    import warnings
    warnings.warn(f"Could not register colormaps due to missing dependencies: {e}", 
                  ImportWarning, stacklevel=2)
except Exception as e:
    import warnings
    warnings.warn(f"Failed to register some colormaps: {e}. "
                  "This may affect colormap functionality.", 
                  RuntimeWarning, stacklevel=2)
