from cgeniepy.array import GriddedData
import numpy as np
import xarray as xr
import pytest


def create_testdata():
    lat = np.linspace(-89.5,89.5,180)
    lon = np.linspace(0,359,360)
    np.random.seed(12349)
    data = np.random.rand(lat.size,lon.size)
    xdata = xr.DataArray(data, coords=[('lat',lat),('lon',lon)],
                         attrs={'long_name':'random data', 'units':'uniteless'})   
    return GriddedData(xdata, attrs=xdata.attrs)

def test_mean():
    data = create_testdata()
    assert data.mean().data.item() == pytest.approx(0.5021667118489245)

def test_std():
    data = create_testdata()
    assert data.std().data.item() == pytest.approx(0.28845163932542145)

def test_std_dim():
    data = create_testdata()
    result = data.std(dim="lat")
    expected = data.data.std(dim="lat")
    xr.testing.assert_allclose(result.data, expected)

def test_sd_alias():
    data = create_testdata()
    assert data.sd().data.item() == pytest.approx(0.28845163932542145)

def test_variance():
    data = create_testdata()
    assert data.variance().data.item() == pytest.approx(0.08320434822952301)

def test_add_and_subtract_align_coordinates():
    data = create_testdata()
    ## the same field stored in a different longitude order, under another name
    rolled = GriddedData(data.data.roll(lon=90, roll_coords=True),
                         attrs={'long_name': 'rolled data', 'units': 'uniteless'})

    xr.testing.assert_allclose((data + rolled).data.sortby('lon'), 2 * data.data)
    xr.testing.assert_allclose((data - rolled).data.sortby('lon'), 0 * data.data)

def test_add_and_subtract_need_the_same_units():
    data = create_testdata()
    other = GriddedData(data.data.copy(), attrs={'units': 'mol/kg'})
    with pytest.raises(ValueError, match="units|mol"):
        data + other
    with pytest.raises(ValueError):
        data - other

def test_power_returns_a_new_object():
    data = create_testdata()
    original = data.data.copy()
    squared = data ** 2

    assert squared is not data
    xr.testing.assert_allclose(data.data, original)
    xr.testing.assert_allclose(squared.data, original ** 2)

def test_median():
    data = create_testdata()
    assert data.median().data.item() == pytest.approx(0.5036269709952814)


def test_min():
    data = create_testdata()
    assert data.min().data.item() ==  pytest.approx(1.2270557283589056e-06)

def test_max():
    data = create_testdata()
    assert data.max().data.item() == pytest.approx(0.9999894976042956)

def test_search_point():
    data = create_testdata()
    ## nemo point lat/lon
    lat = -48.876
    lon = 123.393
    assert data.search_point((lat,lon), ignore_na=True) == pytest.approx(0.25995689209624817)


def create_curvilinear_testdata():
    depth = xr.DataArray([0.0, 10000.0], dims="z_t", attrs={"units": "cm"})
    tlat = xr.DataArray(
        [[0.0, 0.2], [10.0, 10.2]], dims=("nlat", "nlon"), name="TLAT"
    )
    tlong = xr.DataArray(
        [[359.8, 1.0], [0.0, 2.0]], dims=("nlat", "nlon"), name="TLONG"
    )
    values = np.array([
        [[1.0, np.nan], [3.0, 4.0]],
        [[5.0, 6.0], [7.0, 8.0]],
    ])
    array = xr.DataArray(
        values,
        dims=("z_t", "nlat", "nlon"),
        coords={"z_t": depth, "TLAT": tlat, "TLONG": tlong},
    )
    return GriddedData(array)


def test_search_point_curvilinear():
    data = create_curvilinear_testdata()
    assert data.search_point((0.0, 0.0, 0.0)) == pytest.approx(1.0)


def test_search_point_curvilinear_ignore_na_and_depth_units():
    data = create_curvilinear_testdata()
    # The geographically closest surface cell is NaN. At 100 m depth that
    # same cell is valid and closer than distant surface cells.
    assert data.search_point((0.0, 0.2, 1.0), ignore_na=True) == pytest.approx(6.0)


def test_search_point_curvilinear_uses_cf_metadata():
    y = xr.DataArray(
        [[0.0, 0.0], [5.0, 5.0]],
        dims=("j", "i"),
        attrs={"standard_name": "latitude", "units": "degrees_north"},
    )
    x = xr.DataArray(
        [[20.0, 21.0], [20.0, 21.0]],
        dims=("j", "i"),
        attrs={"standard_name": "longitude", "units": "degrees_east"},
    )
    array = xr.DataArray(
        [[1.0, 2.0], [3.0, 4.0]],
        dims=("j", "i"),
        coords={"model_y": y, "model_x": x},
    )
    assert GriddedData(array).search_point((4.9, 20.9)) == pytest.approx(4.0)


def test_search_point_curvilinear_accepts_coordinate_names():
    array = xr.DataArray(
        [[1.0, 2.0], [3.0, 4.0]],
        dims=("row", "column"),
        coords={
            "my_y": (("row", "column"), [[0.0, 0.0], [5.0, 5.0]]),
            "my_x": (("row", "column"), [[20.0, 21.0], [20.0, 21.0]]),
        },
    )
    result = GriddedData(array).search_point(
        (4.9, 20.9), lat_coord="my_y", lon_coord="my_x"
    )
    assert result == pytest.approx(4.0)

def test_sel_modern_basin():
    data = create_testdata()
    assert data.sel_modern_basin(50,norm_lon_method='').mean().data.item() == pytest.approx(0.5019781051132972)
