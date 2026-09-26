from cgeniepy.grid import GridOperation as go
import numpy as np
def test_lon_n2g():
    assert go().lon_n2g(100, -270) == -260


def test_lon_g2n():
    assert go().lon_g2n(-260) == 100


def test_lon_offset_0():
    ## e.g. a grid with par_grid_lon_offset = 0 has longitudes in (0, 360)
    assert go().lon_g2n(355) == -5
    assert go().lon_n2g(-5, grid_lon_offset=0) == 355
    assert go().lon_n2g(5, grid_lon_offset=0) == 5


def test_xr_n2g_offset():
    import xarray as xr
    data = xr.DataArray([1, 2], dims="lon", coords={"lon": [-5, 5]})
    assert list(go().xr_n2g(data, grid_lon_offset=0).lon.values) == [5, 355]


def test_lon_e2n():
    assert go().lon_e2n(350) == -10

def test_lon_n2e():
    assert go().lon_n2e(-10) == 350


def test_geodistance_2d():
    pnt1=(0, 0,0)
    pnt2= np.array([[0, 10, 10]])

    assert go().geo_dis2d(pnt1, pnt2).item() == 1111.9492664455872

def test_geniebin_depth():
    edges = go().get_genie_depth(edge=True)[::-1]
    depths = go().get_genie_depth(edge=False)[::-1]
    assert go().geniebin_depth(0) == depths[0]
    assert go().geniebin_depth(100) == depths[1]
    assert go().geniebin_depth(edges[1]) == depths[1]
    assert go().geniebin_depth(5000) == depths[-1]

def test_checkdimension():
    input = ['lat','lon','time','dpeth'] ## intentional typo
    has_lat, has_lon, has_depth, has_time = go().check_dimension(input)
    assert has_depth == False


def test_geniebin_other_grids():
    assert go().geniebin_lat(10, N=18) in go().get_genie_lat(N=18)
    assert go().geniebin_lon(5, N=18) == 10
    assert go().normbin_lon(5, N=18) == 10
    ## a 16-level grid down to 5500 m
    depths = go().get_genie_depth(N=16, max_depth=5500)[::-1]
    assert go().geniebin_depth(5500, max_depth=5500) == depths[-1]
    assert go().geniebin_depth(100, N=8) == go().get_genie_depth(N=8)[::-1][0]


def test_set_coordinates_3d():
    class Obj: pass
    for index in (['depth', 'lat', 'lon'], ['lat', 'lon', 'depth']):
        obj = Obj()
        go.set_coordinates(obj, index)
        assert (obj.depth, obj.lat, obj.lon) == ('depth', 'lat', 'lon')

    obj = Obj()
    go.set_coordinates(obj, ['time', 'depth', 'lon'])
    assert (obj.time, obj.depth, obj.lon) == ('time', 'depth', 'lon')


def test_dimorder():
    input = ['lat','lon','time','depth']
    depth_order = go().dim_order(input)[1]
    assert input[depth_order] == 'depth'

def test_genie_depth_matches_model():
    import cgeniepy
    model = cgeniepy.sample_model()
    np.testing.assert_allclose(go().get_genie_depth()[::-1], model.grid_mask_3d().data.zt.values)
    np.testing.assert_allclose(go().get_genie_depth(edge=True)[::-1], model.grid_zt_edges().data.values, atol=1e-9)


def test_genie_depth_extra_levels():
    ## 16 levels to 5000 m, plus one level below
    edges = go().get_genie_depth(N=17, edge=True, extra_levels=1)[::-1]
    np.testing.assert_allclose(edges[:17], go().get_genie_depth(edge=True)[::-1], atol=1e-9)
    assert edges[-1] > 5000


def test_genie_lat_equal_degree():
    np.testing.assert_allclose(go().get_genie_lat(N=18, edge=True, equal_area=False), np.arange(-90, 91, 10))
    np.testing.assert_allclose(go().get_genie_lat(N=18, equal_area=False), np.arange(-85, 90, 10))


def test_mask_arctic_med_warns_off_the_modern_grid():
    import warnings
    import pytest
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        go().mask_Arctic_Med(np.ones((36, 36)))
    with pytest.warns(UserWarning, match="modern 36x36"):
        go().mask_Arctic_Med(np.ones((36, 72)))
