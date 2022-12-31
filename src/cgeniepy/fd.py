## A simple functional diversity class, with functional trait as dimension
## (thus not PCoA based) and abundance as weight


import itertools
import numpy as np
from .utils import distance
from .core import GenieArray

class Species:
    "Species class, representing biological species, or functional groups in the biogeochemical modelling"

    def __init__(self, abundance, trait_dict, sub_sp=None):
        "abundance can be a single number or an array"
        self.abundance = abundance
        self.trait = trait_dict
        self.sub_sp = sub_sp

        self.abun_ndim = np.array(self.abundance).ndim
        self.abun_shape = np.array(self.abundance).shape
        self.trait_names = sorted(list(self.trait.keys()))

    def get_trait(self, trait_name: str):
        return self.trait[trait_name]

class Community:

    def __init__(self, sp_list: list):
        self.sp_list = sp_list
        ## Functional entity (group) richness
        self.richnesss = len(sp_list)

        # add into species pool, counting from 1
        for i, sp in enumerate(sp_list):
            setattr(self, f'sp{i+1}', sp)

        self._check_trait()
        self._check_datashape()

        # inherit some attributes
        self.abun_ndim=self.sp_list[0].abun_ndim
        self.abun_shape=self.sp_list[0].abun_shape
        self.trait_names =self.sp_list[0].trait_names

        # handle sub_sp case
        all_sub_sp = [i.sub_sp for i in self.sp_list]

        ## Functional redundancy (average species in each functional group)
        ## i.e., the function number relative to species number
        if None in all_sub_sp:
            self.total_sub_sp = None
        else:
            self.total_sub_sp = sum(all_sub_sp)
            self.redundancy = self.total_sub_sp/self.richnesss

            ## Functional over-redundancy and vulnerability
            for i, sp in enumerate(sp_list):
                # over-redunancy: excessive species' proportion
                over_redundancy = (max(sp.sub_sp, self.redundancy) - self.redundancy)/self.total_sub_sp
                setattr(self.sp_list[i], "over_redundancy", self._init_array(over_redundancy))
                # vulnerability: species proportion in each group
                vulnerability = sp.sub_sp/self.total_sub_sp
                setattr(self.sp_list[i], "vulnerability", self._init_array(vulnerability))

    def _init_array(self, value):
        return np.full(self.abun_shape, value)

    def _check_trait(self):
        "if any two species have different trait, then raise error"
        # check attributes
        for sp1, sp2 in itertools.combinations(self.sp_list, 2):
            if sp1.trait.keys() != sp2.trait.keys():
                raise ValueError("Inconsistent trait names among species!")

    def _check_datashape(self):
        for sp1, sp2 in itertools.combinations(self.sp_list, 2):
            if sp1.abun_shape != sp2.abun_shape:
                raise ValueError("Inconsistent abundance data size among species!")

    def weighted_over_redun(self):
        "Functional over redundancy (only used for same species), unit: redundant species"
        x = self._init_array(0.0)

        for sp in self.sp_list:
            x += sp.abundance * sp.over_redundancy

        g = GenieArray()
        g.array = x
        return g

    def cwm(self, trait: str):
        """community weighted mean, usually weighted in relative abundance
        Geometrically represent the position of community centroid (i,j)
        """

        x = self._init_array(0.0)
        y = self._init_array(0.0)

        for sp in self.sp_list:
            x += sp.get_trait(trait) * sp.abundance
            y += sp.abundance

        # weighted_trait_list / total_abundance_list
        cwm = np.divide(x, y)

        g = GenieArray()
        g.array = cwm
        return g

    def specialisation(self):
        "Functional specialisation, distance bewteen the origin and the centroid"

        origin_tuple = []
        centroid_tuple = []

        for trait in self.trait_names:
            origin_tuple.append(self._init_array(0.0))
            centroid_tuple.append(self.cwm(trait).pure_array())

        fspec = distance(tuple(origin_tuple), tuple(centroid_tuple))

        g = GenieArray()
        g.array = fspec
        return g

    # def fdis(self, *args, **kwargs):
    # "functional dispersion"
    # weight (relative abundance/biomass)
    # weight = self.select_foram("bn").biomass_c()
    # each group's distance to centroid
    # dist =
    # weighted sum

def modern_foram_community(bn_abundance, bs_abundance, sn_abundance, ss_abundance):
    bn = Species(bn_abundance, {"symbiosis": 0, "spine":0}, sub_sp = 17)
    bs = Species(bs_abundance, {"symbiosis": 0, "spine":1}, sub_sp = 2)
    sn = Species(sn_abundance, {"symbiosis": 0.5, "spine":0}, sub_sp = 5)
    ss = Species(ss_abundance, {"symbiosis": 1, "spine":1}, sub_sp = 23)

    foram_community = Community([bn, bs, sn, ss])
    return foram_community    
