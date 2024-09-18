from cgeniepy.skill import ArrComparison
import numpy as np

def create_testdata():
    x = np.linspace(0,100,100)
    y = x    

    ## calculate skill score
    return ArrComparison(x, y)    

def test_mscore():    
    ac = create_testdata()
    assert ac.mscore()==1.0

def test_pearson_r():
    ac = create_testdata()
    diff = ac.pearson_r().item() - 1.0
    assert diff < 1E-8

def test_cos_sim():    
    ac = create_testdata()
    assert ac.cos_similarity()==1.0

def test_rmse():
    ac = create_testdata()
    assert ac.rmse()==0.0
