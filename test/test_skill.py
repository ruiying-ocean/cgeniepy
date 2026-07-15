import matplotlib.pyplot as plt
import numpy as np

from cgeniepy.skill import ArrComparison, TaylorDiagram

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


def test_taylor_diagram_default_colormap():
    comparisons = [
        ArrComparison(np.arange(5), np.arange(5), label="a"),
        ArrComparison(np.arange(5) + 1, np.arange(5), label="b"),
    ]
    diagram = TaylorDiagram(comparisons)
    diagram.setup_ax()
    diagram.plot(add_legend=False)

    assert len(diagram.ax.collections) == len(comparisons)
    plt.close(diagram.fig)
