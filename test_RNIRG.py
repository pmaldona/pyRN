#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 10:21:39 2022

@author: pmaldona
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 11:39:25 2021

@author: pmaldona
"""
#load of the library, it must be in the same work directory
from pyRN import pyRN

# The library contains random network generating functions 
# the most simple random network generator, Nr reactions (>1), Ns species (>1)
# a minimal reaction network is randomly created where each reaction has one reactant and one product and each species
# is used at least once as reactant and once as product, an extra proportion of assignments are carried out randomly
# there are no inflow or outflow reactions, redundant or null reactions may be generated (rn.merge will filter them out)
# dist is a log scaled distribution in the [-1,1] range representing locality
# pr and pp are a log scaled penalization for the repeated use of species as reactants or products
RN = pyRN.rg_g1()

nt=RN.display_RN()
nt.show("RG_0.html")

# New random species can be added to the current reaction network 
# to the existing reactions
RN.rg_extra1()
nt=RN.display_RN()
nt.show("Rg_1.html")

# Adding inflow and otflow reactions
RN.rg_extra_inflow(extra=0.2)
RN.rg_extra_outflow(extra=0.2)
nt.show("RG_2.html")

# Second random generation function which considers the extra_inflow and extra_outflow functions
RN = pyRN.rg_g2()
nt=RN.display_RN()
nt.show("RG_3.html")