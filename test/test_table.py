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