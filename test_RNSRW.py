#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 16:22:45 2022

@author: pmaldona
"""

#load of the library, it must be in the same work directory
from pyRN import pyRN
import numpy as np


# initialization of a randmo generated network:
RN=pyRN.rg_g1(25)    

# Use of the RNSRW module:

# Considering an initialized reaction network an mass action kinetics model
# can be initialized for dynamical simulations proposes be use of the function:
RN.set_model_ma()
# If the function receive an empty argument the initial values are randomly
# initialized from a uniform distribution between [0,1] 
# Then simulation fo 25 time steps can be easily run:
RN.run_model(0,25,100)
# This create a two dataframes, one for the concentration of species and another
# for the reaction rates 
RN.con.plot()
RN.rate.plot()

# The parameters can be changed by means of the following function, 
# where the new concentrations i_sp and kinetic constants of the model rt are 
# entered. The time where the change occurs corresponds to the time at which 
# the simulation was left, or it can be returned by the initial condition init=True.
# At this time only mass action kinetics models can be parameterized.
i_sp=np.random.uniform(size=RN.mr.shape[0])
rt=np.random.uniform(size=RN.mr.shape[1])
RN.param_change_ma(i_sp=i_sp,rt=rt,init=False)
RN.run_model(25,50,100)
RN.con.plot()
RN.rate.plot()



# It is also possible to generate random walks for mass action models.
for i in range(1000):
    RN=pyRN.rg_g1(10)   
    RN.scr_gen_and_pert(w=1,l=10)
