from re import findall
class Chemistry:

    def molar_mass(self, element: str):
        "get molar mass (g/mol) of chemistry elements"
        d = {
            "H": 1.00794,
            "He": 4.002602,
            "Li": 6.941,
            "Be": 9.012182,
            "B": 10.811,
            "C": 12.0107,
            "N": 14.0067,
            "O": 15.9994,
            "F": 18.9984032,
            "Ne": 20.1797,
            "Na": 22.98976928,
            "Mg": 24.305,
            "Al": 26.9815386,
            "Si": 28.0855,
            "P": 30.973762,
            "S": 32.065,
            "Cl": 35.453,
            "Ar": 39.948,
            "K": 39.0983,
            "Ca": 40.078,
            "Sc": 44.955912,
            "Ti": 47.867,
            "V": 50.9415,
            "Cr": 51.9961,
            "Mn": 54.938045,
            "Fe": 55.845,
            "Co": 58.933195,
            "Ni": 58.6934,
            "Cu": 63.546,
            "Zn": 65.409,
            "Ga": 69.723,
            "Ge": 72.64,
            "As": 74.9216,
            "Se": 78.96,
            "Br": 79.904,
            "Kr": 83.798,
            "Rb": 85.4678,
            "Sr": 87.62,
            "Y": 88.90585,
            "Zr": 91.224,
            "Nb": 92.90638,
            "Mo": 95.94,
            "Tc": 98.9063,
            "Ru": 101.07,
            "Rh": 102.9055,
            "Pd": 106.42,
            "Ag": 107.8682,
            "Cd": 112.411,
            "In": 114.818,
            "Sn": 118.71,
            "Sb": 121.760,
            "Te": 127.6,
            "I": 126.90447,
            "Xe": 131.293,
            "Cs": 132.9054519,
            "Ba": 137.327,
            "La": 138.90547,
            "Ce": 140.116,
            "Pr": 140.90465,
            "Nd": 144.242,
            "Pm": 146.9151,
            "Sm": 150.36,
            "Eu": 151.964,
            "Gd": 157.25,
            "Tb": 158.92535,
            "Dy": 162.5,
            "Ho": 164.93032,
            "Er": 167.259,
            "Tm": 168.93421,
            "Yb": 173.04,
            "Lu": 174.967,
            "Hf": 178.49,
            "Ta": 180.9479,
            "W": 183.84,
            "Re": 186.207,
            "Os": 190.23,
            "Ir": 192.217,
            "Pt": 195.084,
            "Au": 196.966569,
            "Hg": 200.59,
            "Tl": 204.3833,
            "Pb": 207.2,
            "Bi": 208.9804,
            "Po": 208.9824,
            "At": 209.9871,
            "Rn": 222.0176,
            "Fr": 223.0197,
            "Ra": 226.0254,
            "Ac": 227.0278,
            "Th": 232.03806,
            "Pa": 231.03588,
            "U": 238.02891,
            "Np": 237.0482,
            "Pu": 244.0642,
            "Am": 243.0614,
            "Cm": 247.0703,
            "Bk": 247.0703,
            "Cf": 251.0796,
            "Es": 252.0829,
            "Fm": 257.0951,
            "Md": 258.0951,
            "No": 259.1009,
            "Lr": 262,
            "Rf": 267,
            "Db": 268,
            "Sg": 271,
            "Bh": 270,
            "Hs": 269,
            "Mt": 278,
            "Ds": 281,
            "Rg": 281,
            "Cn": 285,
            "Nh": 284,
            "Fl": 289,
            "Mc": 289,
            "Lv": 292,
            "Ts": 294,
            "Og": 294,
            "": 0,
        }

        return d[element]


    def _formula_parser(self, formula: str):
        """
        a chemical formular parser

        :param formula: formula string
        :returns: dictionary with element (key) and number (value)
        """
        s = findall("([A-Z][a-z]?)([0-9]*)", formula)
        d = dict(s)
        for k, v in d.items():
            if v == "":
                d[k] = 1
            else:
                d[k] = int(v)
        return d


    def molecular_weight(self, formula: str):
        """
        Calculate molecular weight using a chemical formula

        :param formula: chemical formula
        :returns: molecular weight

        ---------
        Example
        ----------
        >>> molecular_weight("C6H12O6")
        >>> 180.15588
        """
        "calculate molecular weight"
        d = self._formula_parser(formula)
        mw = sum([self.molar_mass(k) * v for k, v in d.items()])
        return mw


    def pure_unit(self, unit: str):
        """
        Get a pure unit string by removing elements in string (e.g., mmol C d-1)

        :param unit: unit string with elements
        :returns: pure unit string

        ---------
        Example
        ----------
        >>> pure_unit("mmol C d-1")
        >>> "mmol d-1"
        """
        # an element pool
        l = ["C", "N", "P", "O", "Fe", "Si", "Ca", "S"]
        # add leading whitespace
        l = [" " + i for i in l]

        # remove element in unit string
        for element in l:
            if element in unit:
                unit = unit.replace(element, "")

        return unit


    def format_base_unit(self, input_str, latex_format=True):
        """
        Convert the base unit (i.e., no whitespace) to the desired format
        
        base unit means it is not a compound unit (e.g., m2, m-2, m/s, m s-1),
        which should be handled by the function `format_unit`
                
        :param input_str: a base unit string 
        :return: a formatted string
        
        ----------
        Example:
        ----------
        
        >>> convert_a_unit('m2')
        '$m^2$'
        >>> convert_a_unit('m-2')
        '$m^{-2}$'
        """

        if input_str == 'n/a':
            return None

        replacements = {
            '^': '',
            '**': '^',
            'o/oo': '‰',
            'degrees': '°',
            'degC': '°C',
            '° C': '°C',
        }
        
        ## replace some special characters
        for k, v in replacements.items():
            if k in input_str:
                input_str = input_str.replace(k, v)

        ## when exponent is negative  (e.g., m-2, m-3)
        if '-' in input_str:
            parts = input_str.split('-')
            ## two parts: base and exponent
            base = parts[0]
            exponent = parts[1]
            
            if len(parts) == 2 and exponent.isdigit():
                exponent = int(exponent)
                if exponent != 0:
                    # Convert to the desired format using string formatting
                    # 3 layers of curly braces are needed to print a single curly brace
                    output_str = f"{parts[0]}$^{{-{exponent}}}$"
                    ## if not use LaTeX format, then remove dollar sign
                    if not latex_format: output_str = output_str.replace('$', '')                
                    return output_str
                
        ## when exponent is positive (e.g., m2, m3)
        ## also slice into two parts: base and exponent
        if len(input_str) > 1:
            base = input_str[:-1]
            exponent = input_str[-1]
            
            if exponent.isdigit():
                exponent = int(exponent)
                if exponent != 0:
                    # Convert to the desired format using string formatting
                    # 3 layers of curly braces are needed to print a single curly brace
                    output_str = f"{base}$^{{{exponent}}}$"
                    ## if not use LaTeX format, then remove dollar sign
                    if not latex_format: output_str = output_str.replace('$', '')
                    return output_str
        
        return input_str

    def format_unit(self, input, *args, **kwargs):
        """
        convert a unit in a string seprated by whitespace
        to the desired format, which can be used in plot labels and 
        facilitate the unit calculation
        
        :param input: a unit string 
        :return: a  formatted string
        
        ----------
        Example:
        ----------
        
        >>> format_unit('mol kg-1')
        'mol kg$^{-1}$'
        >>> format_unit('mmol m-2 d-1')
        'mmol m$^{-2}$ d$^{-1}$'
        """
        if input == 'n/a':
            return None

        if '{' in input or '}' in input:
            return input

        if input == '° C':
            return '°C'
        
        unit_list = input.split(" ")
        unit_list = [self.format_base_unit(i, *args, **kwargs) for i in unit_list]

        clean_string = ' '.join(unit_list)
        return clean_string



class ModelCarbon:
    
    def __init__(self, model):
        self.model = model
    
    def get_NPP(self):
        "netprimary production (NPP) in Pg C yr-1"
        npp = self.model.get_var("eco2D_Uptake_Fluxes_C").isel(time=-1)
        grid_vol_sur = self.model.grid_volume().isel(time=-1).isel(zt=0) * 1E-3
        npp = npp * grid_vol_sur
        npp = npp * 12 * 1e-15 * 365
        return npp

    def get_poc_export(self, depth=None):
        "particulate organic carbon (POC) export at 80 m in Pg C yr-1"
        poc_export_80 = self.model.get_var("bio_fexport_POC").isel(time=-1) * 12.01 / 1e15

        if not depth:
            return poc_export_80
        else:
            # if depth is specified, return POC export at the specified depth
            model_area = self.model.grid_area()        
            fpart_poc = self.model.get_var("bio_fpart_POC").isel(time=-1) * 12.01 / 1e15 * model_area
        
            poc_export_deep = fpart_poc.sel(zt=depth, method='nearest')
            return poc_export_deep

    def get_eratio(self):
        "export ratio (ER), calculated as POC export at 80 m divided by NPP"
        poc_export_80 = self.model.get_var("bio_fexport_POC").isel(time=-1) * 12.01 / 1e15    
        npp = self.get_NPP()
        eratio = poc_export_80 / npp
        return eratio

    def get_Teff(self, zt1=80, zt2=1000):
        """
        transfer efficiency (Teff), default calculated as POC export at 1000 m divided by POC export at 80 m
        """
        model_area = self.model.grid_area()        
        fpart_poc = self.model.get_var("bio_fpart_POC").isel(time=-1) * 12.01 / 1e15 * model_area
        
        poc_export_shallow = fpart_poc.sel(zt=zt1, method='nearest')
        poc_export_deep = fpart_poc.sel(zt=zt2, method='nearest')

        transfer_efficiency = poc_export_deep / poc_export_shallow
        return transfer_efficiency
