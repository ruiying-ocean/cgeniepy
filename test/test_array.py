from cgeniepy.array import GriddedData
import numpy as np
import xarray as xr

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
    assert data.mean().data.item() == 0.5021667118489245

def test_sd():
    data = create_testdata()
    assert data.sd().data.item() == 0.28845163932542145

def test_variance():
    data = create_testdata()
    assert data.variance().data.item() == 0.08320434822952301

def test_median():
    data = create_testdata()
    assert data.median().data.item() == 0.5036269709952814


def test_min():
    data = create_testdata()
    assert data.min().data.item() == 1.2270557283589056e-06    

def test_max():
    data = create_testdata()
    assert data.max().data.item() == 0.9999894976042956

def test_search_point():
    data = create_testdata()
    ## nemo point lat/lon
    lat = -48.876
    lon = 123.393
    assert data.search_point((lat,lon), ignore_na=True) == 0.25995689209624817

def test_sel_modern_basin():
    data = create_testdata()
    assert data.sel_modern_basin(50,norm_lon_method='').mean().data.item() == 0.5019781051132972
