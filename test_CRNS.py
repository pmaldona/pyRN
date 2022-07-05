#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 11:39:25 2021

@author: pmaldona
"""
#load of the library, it must be in the same work directory
from CRNS import CRNS
import time
import networkx as nx
import numpy as np
from bitarray import bitarray as bt

# loading form a text file please refer to rn_test.txt to see example
file="networks/rn_test.txt"
RN = CRNS.from_txt(file)

# file="networks/PW000035.sbml"
# file="networks/Farm_sbml.xml"
# RN = CRNS.from_sbml(file,False)

# Basic Variables
print(RN.mp)
print(RN.mr)
print(RN.sp)
print(RN.sp_n)


start = time.time()
# generation of the basics sets
b_steps=RN.gen_basics()
end = time.time()
b_time=end-start

# Plot histogram of species and reactions related to basics and partitions 
RN.plot_basic_r_presence([0,2,3])
RN.plot_basic_sp_presence(RN.sp[[0,2,3]])

start = time.time()
# generation of the strucutres close structre
syn_steps=RN.gen_syn_str()
end = time.time()
syn_time=end-start

# generation of all proto-synergies
RN.gen_mgen()
RN.all_syn()
# finding all synergies for an especific organization
RN.syn_sets(RN.syn_org[4])


start = time.time()
# Generation of the semi-self-maintained connected structure
ssm_steps=RN.gen_ssm_str()
end = time.time()
ssm_time=end-start


start = time.time()
# Generation of the dynamical connected structure
dyn_steps=RN.gen_dyn_str()
end = time.time()
dyn_time=end-start


start = time.time()
# Hasse of overproduced species
dyn_dcom=[]
op_h_steps=[]
for i in RN.dyn_org:
    if len(i)>1:
        G , op_h_st = RN.over_hasse(i)
        dyn_dcom+=list(nx.get_node_attributes(G,'dcom').values())
        op_h_steps.append(op_h_st)
    else:
        dyn_dcom.append(RN.dcom(i, i, RN.sp2r(i)))

end = time.time()
op_h_time=end-start

steps=[b_steps,syn_steps,ssm_steps,dyn_steps,op_h_steps]
times=[b_time,syn_time,ssm_time,dyn_time,op_h_time]

# Network dynamic test 

# loading mass action model
RN.ma_model()
# simulation fo 25 time steps
print(RN.model.simulate(0,25,100))
# plot of the dynamic
RN.model.plot()

# Loading sbml file
# RN=CRNS.from_sbml("../networks/BIOMD0000000013.xml")
RN=CRNS.from_sbml("../COT/networks/ReacNet/raw/BIOMD0000000013.xml")
# loading mass action model
RN.ma_model(i_sp=np.ones(RN.mp.shape[0]),rt=np.ones(RN.mp.shape[1]))
# simulation fo 25 time steps
print(RN.model.simulate(0,25,100))
# plot of the dynamic
RN.model.plot()

#loading sbml model
RN.sbml_model()
# simulation fo 25 time steps
print(RN.model.simulate(0,25,100))
# plot of the dynamic
RN.model.plot()
