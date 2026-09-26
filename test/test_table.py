from cgeniepy.table import ScatterData
import pandas as pd

def create_testdata():
    lat = -48.876
    lon = 123.393
    df = pd.DataFrame({'lat': [lat], 'lon': [lon]})
    return ScatterData(df)

def test_init():
    data = create_testdata()
    ## note the order is put in reverse intentionally
    data.set_index(['lon','lat']) 
    assert data.lat == 'lat'

def test_detectbasin():
    data = create_testdata()
    data.set_index(['lat','lon'])
    basin_value = data.detect_basin()['basin'].values.item()
    assert basin_value =='Indian Ocean'

def test_url():
    url = "https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/epica_domec/edc-co2-2008-bern-noaa.txt"
    
    test_data= ScatterData(url, comment='#', delimiter='\t')
    assert test_data.data['CO2'][0] == 257.8

def test_to_geniebin_model_grid():
    import cgeniepy
    model = cgeniepy.sample_model()
    zt = model.grid_mask_3d().data.zt.values
    lat = model.grid_mask().data.lat.values
    df = pd.DataFrame({'depth': [100., 5000., 6000.], 'lat': [0.5, -90., 10.],
                       'lon': [179., -180., 10.], 'v': [1., 2., 3.]})
    data = ScatterData(df)
    data.set_index(['depth', 'lat', 'lon'])
    binned = data.to_geniebin('v', model=model).reset_index()
    ## 6000 m is below the ocean floor of the model
    assert len(binned) == 2
    assert set(binned['depth']) == {zt[1], zt[-1]}
    assert set(binned['lat']) <= set(lat)
    assert set(binned['lon']) == {175.}
