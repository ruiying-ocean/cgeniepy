import importlib
import warnings

import cgeniepy
import matplotlib.colors

from cgeniepy.plot import CommunityPalette


def test_txt_palette():
    cmap = CommunityPalette().get_palette('ODV')
    assert matplotlib.colors.to_hex(cmap.colors[0]) == '#feb483'


def test_xml_palette():
    cmap = CommunityPalette().get_palette('Section')
    assert matplotlib.colors.to_hex(cmap.colors[0]) == '#c2ccb8'


def test_alt_init():
    hex_codes = CommunityPalette('my_rainbow').to_hex()
    assert hex_codes[0] == '#320064'


def test_colormap_registration_is_idempotent():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        importlib.reload(cgeniepy)


def test_phanda_palette():
    """the coolors.co teal-to-red palette used for the climate states in Judd et al. 2024 (Science)"""
    hex_codes = CommunityPalette('PhanDA').to_hex()
    assert hex_codes[0] == '#001219' and hex_codes[-1] == '#9b2226'
    assert 'PhanDA' in matplotlib.colormaps


def test_tlylrd_palette():
    """PhanDA's hues rebuilt as a diverging map: lightness linear on each arm and mirrored about the sand centre"""
    hex_codes = CommunityPalette('TlYlRd').to_hex(unique=False)
    assert len(hex_codes) == 256
    assert hex_codes[0] == '#006173' and hex_codes[-1] == '#9a1e1e'
    ## symmetric: the lightest colour sits at the midpoint
    brightest = max(range(256), key=lambda i: sum(matplotlib.colors.to_rgb(hex_codes[i])))
    assert brightest in (127, 128)
    assert 'TlYlRd' in matplotlib.colormaps
