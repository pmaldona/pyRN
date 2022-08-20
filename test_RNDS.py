#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:12:02 2022

@author: pmaldona
"""
#load of the library, it must be in the same work directory
from pyRN import pyRN
import networkx as nx

# initialization of a randmo generated network:
RN=pyRN.rg_g2()    

# Use of the RNDS module:

# For a given set of species and set of reactions, the the optimal number of 
# overproduced species can be obtained as follow:
op_sp=RN.over_prod(RN.syn_org[0], RN.sp2r(RN.syn_org[0]))
print(RN.bt_to_sp(op_sp))
# Or the decomposition it self, considering an imput of overproduced species
dcom=RN.dcom(op_sp, RN.syn_org[0], RN.sp2r(RN.syn_org[0]))
print(dcom)
# The value of the dcom vector is -1 it corresponds to an 
# overproducible species, if it is -2 to a catalytic species and if it 
# is 0 the species is not present. The integer values indicate belonging 
# to the current fragile cycle.

# Also all possible decomposition for a given set can be obtained, this is done
# by use of the following function:
G , op_h_st = RN.over_hasse(RN.syn_org[0])    

# G corresponds to a networkx graph that is the Hasse diagram of the overproducible
# species, each node has as attribute the decomposition, op_h_sp, are the number 
# of steps to create the Hasse diagram. 
all_dcom = list(nx.get_node_attributes(G,'dcom').values())
print(all_dcom)
# For more details of use of this module please refer to the "./pyRN/RNDS.py" 
# file.