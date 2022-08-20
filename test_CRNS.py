#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:11:40 2022

@author: pmaldona
"""
from pyRN import pyRN
import time

# Use of the CRNS module:
# Initi of the example network
file="networks/rn_test.txt"
RN = pyRN.from_txt(file)

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