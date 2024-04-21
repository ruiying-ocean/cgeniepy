import cgeniepy.chem as chem

def test_molcularweight():
    assert chem.Chemistry().molecular_weight('CO2') == 44.0095


def test_formatunit():
    input = 'mmol C m-2 yr-1'
    assert chem.Chemistry().format_unit(input) == 'mmol C m$^{-2}$ yr$^{-1}$'