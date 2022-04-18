#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 11:39:25 2021

@author: pmaldona
"""
#load of the library, it must be in the same work directory
from CRNS import CRNS
from bitarray import frozenbitarray
import networkx as nx
from bitarray import bitarray as bt

# loading form a text file please refer to rn_test.txt to see example
file="networks/rn_test.txt"
RN = CRNS.form_txt(file)

# file="networks/PW000035.sbml"
# RN = CRNS.from_sbml(file,False)

print(RN.mp)
print(RN.mr)
print(RN.sp)
print(RN.sp_n)

# # Veryfing if network is RN-mantainig
print(RN.is_sm(RN.sp))


# # Generation of reaction equivalence clases and basic sets
RN.gen_basics()
print(RN.sp_b)
print(RN.conn)
print(RN.dyn_conn)

# # Generation of the synergic structure
RN.gen_syn_str()
print(RN.syn_str.nodes())
print(RN.syn_str.edges())
print(RN.ssms)

selected_edges = [(u,v,e) for u,v,e in RN.syn_str.edges(data=True) if e['syn'] == True]

# print(RN.syn_str.nodes[frozenbitarray('000101')])
for u,v,e in RN.syn_str.edges(data=True):
    if e['syn']==True:
        print ([u,v,RN.syn_str[u][v]])

# Generation of the synergic structure
RN.gen_ssm_str()
print(RN.ssm_str.nodes())
print(RN.ssm_str.edges())
print(RN.ssms)

# # generation of minimal generators
RN.gen_mgen()
print(len(RN.mgen))

# Generation of all proto synergies
RN.all_syn()  
print(len(RN.syn))
# print(RN.syn_p)


pr=bt(RN.mp.shape[1])
pr.setall(1)
sp=bt(RN.mp.shape[0])
sp.setall(1)

# generating the decomposition of a reaction network
dcom=RN.op_dcom(sp, pr)
print(dcom)