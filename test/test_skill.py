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
    assert ac.pearson_r()==0.9999999999999999

def test_cos_sim():    
    ac = create_testdata()
    assert ac.cos_similarity()==1.0

def test_rmse():
    ac = create_testdata()
    assert ac.rmse()==0.0