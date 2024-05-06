from cgeniepy.grid import GridOperation as go
import numpy as np
def test_lon_n2g():
    assert go().lon_n2g(100, -270) == -260


def test_lon_g2n():
    assert go().lon_g2n(-260) == 100


def test_lon_e2n():
    assert go().lon_e2n(350) == -10

def test_lon_n2e():
    assert go().lon_n2e(-10) == 350


def test_geodistance_2d():
    pnt1=(0, 0,0)
    pnt2= np.array([[0, 10, 10]])

    assert go().geo_dis2d(pnt1, pnt2).item() == 1111.9492664455872

def test_checkdimension():
    input = ['lat','lon','time','dpeth'] ## intentional typo
    has_lat, has_lon, has_depth, has_time = go().check_dimension(input)
    assert has_depth == False


def test_dimorder():
    input = ['lat','lon','time','depth']
    depth_order = go().dim_order(input)[1]
    assert input[depth_order] == 'depth'