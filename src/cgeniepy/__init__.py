# from pint import UnitRegistry
# import pathlib


# ureg = UnitRegistry()
# Q_ = ureg.Quantity
# file_path = pathlib.Path(__file__).parent.parent / "data/context.txt"
# ureg.load_definitions(file_path)


from importlib.resources import files
from cgeniepy.model import GenieModel
from cgeniepy.ecology import EcoModel

def sample_model(model_type='GenieModel', *args, **kwargs):
    file_path=str(files('data').joinpath('sample_model'))
    
    if model_type == 'GenieModel':
        model = GenieModel(file_path,*args, **kwargs)
    elif model_type == 'EcoModel':
        model = EcoModel(file_path,*args, **kwargs)
    else:
        raise ValueError('model_type must be either GenieModel or EcoModel')    

    return model
    
