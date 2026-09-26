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

def test_crmse():
    assert create_testdata().crmse() == 0.0

    model = np.array([1.0, 3.0, 2.0, 5.0, np.nan])
    obs = np.array([2.0, 2.0, 4.0, 3.0, 1.0])
    ac = ArrComparison(model, obs)
    m, o = model[:4], obs[:4]
    sm, so, r = m.std(), o.std(), np.corrcoef(m, o)[0, 1]
    # Taylor (2001): crmse**2 = sm**2 + so**2 - 2*sm*so*r
    assert np.isclose(ac.crmse(), np.sqrt(sm**2 + so**2 - 2 * sm * so * r))
    # a constant offset leaves the centred error unchanged
    assert np.isclose(ArrComparison(model + 10, obs).crmse(), ac.crmse())


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
