from cgeniepy.plot import CommunityPalette
import matplotlib.colors


def test_txt_palette():
    cmap = CommunityPalette().get_palette('ODV')
    assert matplotlib.colors.to_hex(cmap.colors[0]) == '#feb483'


def test_xml_palette():
    cmap = CommunityPalette().get_palette('Section')
    assert matplotlib.colors.to_hex(cmap.colors[0]) == '#c2ccb8'


def test_alt_init():
    hex_codes = CommunityPalette('my_rainbow').to_hex()
    assert hex_codes[0] == '#320064'
