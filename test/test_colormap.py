from cgeniepy.plot import CommunityPalette

def test_get_palette():

    cmap = CommunityPalette().get_palette('ODV')
    assert cmap.colors[0] == '#FEB483'


