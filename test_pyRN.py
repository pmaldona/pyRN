#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 11:39:25 2021

@author: pmaldona
"""
#load of the library, it must be in the same work directory
from pyRN import pyRN
import time
import networkx as nx
import numpy as np

# loading form a text file please refer to rn_test.txt to see example,
# textfiles preserve antimony general strucutre in terms of reaction and 
# reaction arrows.
file="networks/rn_test.txt"
RN = pyRN.from_txt(file)

# Alternative sbml files can be loaded. 
# file="networks/PW000035.sbml"
# RN = pyRN.from_sbml(file,False)

# Basic Variables can easily obtain
print("Vector of species")
print(RN.sp)
print("Vector of species names, in case of an initialization form an smbl file result can be different")
print(RN.sp_n)
print("Display species")
print(RN.sin_print_sp())
print("Display reactions")
print(RN.sin_print_r())
print("Reactive stochiometric matrix")
print(RN.mr)
print("Productive stoichimetric matrix")
print(RN.mp)
print("Stoichiometric Matrix")
print(RN.mp-RN.mr)
RN.plot_S()
nt=RN.display_RN()
nt.show("RN.html")


# Use of the CRNS module:

start = time.time() 
# Generation of the basics sets, this function generates all basic set, 
# it creates list of bitarrays for different variables related to the 
# basics sets. Please refer to ./pyRN/CRNS.py module for details.
b_steps=RN.gen_basics()
end = time.time()
b_time=end-start

# The function returns the number of steps (b_steps) in generating all basics
# set, b_time is real the time to compute.

# The generated variables form the function corresponds to bitarrays, for a proper 
# visualization the function RN.bt_to_sp(). For example, the species of 
# the first basics stored as RN.sp_b[0] variable. This species can be seen as 
# proper species by use the following function:
print("example of a basic set")
print(RN.bt_to_sp(RN.sp_b[0]))

print("We can also obtain the triggering reactions for this set of species")
print(RN.bt_ind(RN.sp2r(RN.sp_b[0])))

# There are also functions that plot histogram of species and reactions 
# related to basics and partitions in which they are contained 
RN.plot_basic_sp_presence(RN.sp[[0,2]])
RN.plot_basic_r_presence([0,2])


start = time.time()
# Generation of the synergistic structure, this function create a networkx 
# multigraph (RN.syn_str) and two list of bitarrays corresponding to the 
# semi-sef-maintained sets (RN.syn_ssms) and organizations (RN.syn_org)
# Please refer to ./pyRN/CRNS.py module for details.
syn_steps=RN.gen_syn_str()
end = time.time()
syn_time=end-start

# As same as b_steps, syn_steps are the number of steps in generating the 
# synergistic structure, ans syn_time the time to compute.

# the synergistic edges can be obtain by searching by the property
syn_edges = [(u,v) for u,v,e in RN.syn_str.edges(data=True) if e['syn']]
print("synergetic edges of the RN")
print(syn_edges)

# Also the synergetic structure can be display
nt=RN.display_str(RN.syn_str)
nt.show("CStr.html")

# There is also two (gen_ssm_str and gen_syn_str) function that generate similar 
# structures. They don't give the complete synergistic structure, it's a reduce 
# structure to aim a quicker calculation to obtain the organizations and the 
# semi-self-maintained set, please refer to "./pyRN/CRNS.py" module for more
# details

# Generation of all proto-synergies, to obtain all proto-synergies, in first instance
# it is necessary to generate the minimal generator for each partition.
RN.gen_mgen()
print("minimal generators of the network")
print(RN.mgen)
# After this all proto-synergies can be generated.   
RN.all_syn()
print("All synergies of the network")
print(RN.syn)
print(RN.syn_p)

# Also the proto-synergies can be display
nt=RN.display_pr_syn()
nt.show("pr_syn.html")

print("new part of code, until here")
# Both gen_mgen and all_syn, generates member lists of bitarrays related to the 
# proto-synergies, for more details please refer to "./pyRN/CRNS.py module .

# There is also an function (RN.syn_sets) that for a given set, calculates all possibles 
# synergies that can be done whit it. This function can only be run after generating all
# proto-synergies, in the case considering the second organization this corresponds to:

    
print(RN.syn_sets(RN.syn_org[1]))

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


# Use of the RNSRW module:

# Considering an initialized reaction network an mass action kinetics model
# can be initialized for dynamical simulations proposes be use of the function:
RN.set_model()
# If the function receive an empty argument the initial values are randomly
# initialized from a uniform distribution between [0,1] 
# Then simulation fo 25 time steps can be easily run:
RN.run_model(0,25,100)
# This create a two dataframes, one for the concentration of species and another
# for the reaction rates 
print(RN.con)
print(RN.rate)

# Sbml files can be loaded for inlcuding specific dynamics definitions. Frist
# as rection network object
# RN=pyRN.from_sbml("../COT/networks/ReacNet/raw/BIOMD0000000013.xml")
# And then as model
RN.set_sbml_model()
#  Afterwards the simulation can be run.
RN.run_model(0,25,100)
print(RN.con)
print(RN.rate)

# Other option is to initialize an MKA model from an sbml file whit initial concentrations
# and rates.
RN.set_model(i_sp=np.ones(RN.mp.shape[0]),rt=np.ones(RN.mp.shape[1]))
RN.run_model(0,25,100)
print(RN.con)
print(RN.rate)

# For more details of use of this module please refer to the "./pyRN/RNSRW.py" 
# file.

