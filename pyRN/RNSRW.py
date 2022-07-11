#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:48:31 2022

@author: pmaldona

Reaction Network Simulator and Random Walk Class
"""
from .RNIRG import RNIRG
import numpy as np
import tellurium as te
import numpy as np
import pandas as pd



class RNSRW(RNIRG):
    
    # Funtions that create a mass action dinamics telurrium model (CRNS.model) of 
    # reaction network. It recives as input a vector of the initial concentration 
    # of species i_sp, the reactive constant vector rt and the concentration 
    # thearshold where a species is condidered not present. if the variables 
    # i_sp and rt are not given, they are randomly initialized by a uniform 
    # distribution between 0 and 1
    def ma_model(self, i_sp=np.array([]) ,rt=np.array([]) ,th=0.1):
        
        # Creating the model
        self.model=te.loada("")
        self.model.addCompartment("C", 1,True)
        
        # Creating the random initial concetration if it's out of condition
        if (i_sp.size==0): 
            i_sp=np.random.rand(self.mp.shape[0])
        elif(len(i_sp) !=self.mp.shape[0]):
            i_sp=np.random.rand(self.mp.shape[0])
        
        # Creating the random initial rates if it's out of condition
        if (rt.size==0):  
            rt=np.random.rand(self.mp.shape[1])
        elif (len(rt) !=self.mp.shape[1]):
           rt=np.random.rand(self.mp.shape[1])
           
        # Adding species to the model and treshold concetration
        for i in range(len(self.sp)):
            self.model.addSpecies(self.sp[i], compartment="C", initConcentration=i_sp[i],forceRegenerate=True)
            self.model.addEvent("e"+str(i),False,self.sp[i]+"<"+str(th),True)
            self.model.addEventAssignment("e"+str(i),self.sp[i],"0",True)
        
        # Adding reactions and reactions rate constants
        for i in range(self.mp.shape[1]):
            self.model.addParameter("k"+str(i),rt[i],True)
            
            prod=[]
            for j in np.where(self.mp.iloc[:,i])[0]:
                if self.mp.iloc[j,i].astype(int) > 1:
                    prod.append(str(self.mp.iloc[j,i].astype(int))+self.mp.index[j])
                else:
                    prod.append(self.mp.index[j])
            reac=[]
            rate="k"+str(i)
            for j in np.where(self.mr.iloc[:,i])[0]:
                if self.mr.iloc[j,i].astype(int) > 1:
                    reac.append(str(self.mr.iloc[j,i].astype(int))+self.mr.index[j])
                else:
                    reac.append(self.mr.index[j])
                    
                if self.mr.iloc[j,i]>1:
                    for k in range(1,self.mr.iloc[j,i].astype(int)+1):
                        rate+=" * "+self.mr.index[j]
                else:
                    rate+=" * "+self.mr.index[j]
   
            self.model.addReaction("r"+str(i), reac, prod, rate)
                    
    # Function that load dynamical model directly form the sbml file     
    def sbml_model(self):
        if not self.sbml:
            print("network do not correspond to an sbml file")
            return
        self.model = te.loadSBMLModel(self.fname)
   
    # Function that simulates the dynamics of the proposed model. It takes as 
    # input the initial time ti and the final time tf, and the number of steps 
    # to simulate by means of the variable steps. The function generates the 
    # dataframes CRNS.with CRNS.rate, which correspond to concentrations and 
    # processes of the simulation respectively.
    def simulate(self,ti=0,tf=50,steps=100): 
        
        if ti <= self.model.getCurrentTime():
            t_step=(tf-ti)/steps
            con=[[ti]+self.model.getFloatingSpeciesConcentrations().tolist()]
            rate=[[ti]+self.model.getReactionRates().tolist()]
            for i in range(steps):
                t=t_step*i
                self.model.oneStep(t,t_step)
                con.append([t+t_step]+self.model.getFloatingSpeciesConcentrations().tolist())
                rate.append([t+t_step]+self.model.getReactionRates().tolist())
            
            con=pd.DataFrame(con,columns=['time']+self.model.getFloatingSpeciesConcentrationIds())
            con=con.set_index('time')
            
            rate=pd.DataFrame(rate,columns=['time']+self.model.getReactionIds())
            rate=rate.set_index('time')
            
            self.con=con
            self.rate=rate
        else:
            print("initial time not simulated yet")
    
    
    # Function that modify parameters of a given model. It recives as input a 
    # vector of the initial concentration of species i_sp, the reactive 
    # constant vector rt and the concentration 
    
    def param_model(self, i_sp=np.array([]) ,rt=np.array([])):
        # Creating the random initial concetration if it's out of condition
        if (i_sp.size==0): 
            i_sp=np.random.rand(self.mp.shape[0])
        elif(len(i_sp) !=self.mp.shape[0]):
            i_sp=np.random.rand(self.mp.shape[0])
        
        # Creating the random initial concetration if it's out of condition
        if (rt.size==0):  
            rt=np.random.rand(self.mp.shape[1])
        elif (len(rt) !=self.mp.shape[1]):
            rt=np.random.rand(self.mp.shape[1])
         
        if self.sbml:
            bc_sp=np.where(list(map(lambda x: x in self.model.getBoundarySpeciesIds(), self.mp.index.tolist())))[0]
            f_sp=np.where(list(map(lambda x: x in self.model.getFloatingSpeciesIds(), self.mp.index.tolist())))[0]
