import cgeniepy

def test_model_getvar():
    model = cgeniepy.sample_model()
    assert model.get_var("ocn_sur_temp").mean().data.values.item() == 18.055744171142578
    
def test_eco_pft():
    model = cgeniepy.sample_model(model_type='EcoModel',gemflag=['ecogem'])
    ## 1st PFT
    data = model.get_pft(1).isel(time=-1).mean().data.values.item()
    assert data  == 0.12259095907211304

def test_eco_multipft():
    model = cgeniepy.sample_model(model_type='EcoModel',gemflag=['ecogem'])
    ## 1st and 2nd PFT
    data = model.get_pft([1,2]).isel(time=-1).sum(dim='variable').mean().data.values.item()
    assert data  == 0.16739805042743683
