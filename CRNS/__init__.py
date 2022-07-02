#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 11:39:25 2021

@author: Pedro Maldonado

# Closed Reaction Network Structure Library

"""
import numpy as np
import pandas as pd
import re
from bs4 import BeautifulSoup as bs
import copy
from bitarray import bitarray as bt
from bitarray import frozenbitarray as fbt
from scipy.optimize import linprog
import networkx as nx
import pypoman as ph
from itertools import chain, combinations
from typing import List, Any, Iterable
import random as rm
import tellurium as te


#Class definition: It consist
class CRNS:
    # the inizalization of the obnject generates an stoichmetric matrix for the 
    # reactants (RN.mr) as for the products (RN.mp). The default inizalizaction
    #  is considered form a text file which represents the reaction network as 
    #  writen reaction similar to the antimony format. Also reaction networks 
    # can be inizalizated from antimony fomrtas, smbl and well as crn   

    def __init__(self):
            self.sp=None
            self.sp_n=None
            self.mr=None
            self.mp=None
            self.prod=None
            self.reac=None

    # Reaction Network inizalization form text file similar to antimony format, 
    # see example "rn_test.txt". calling the object can be done by RN(file_name). 
    # Example rn = RN.form_txt("rn_text.txt")
    @classmethod
    def from_txt(cls,file):
        try:
            # list fitering
            rn = pd.read_csv(file,header=None)
            rn[0] = rn[0].str.replace("[ \t]","",regex=True)        
            rn[0] = rn[0].str.replace("^.*:","",regex=True)        
            rn[0] = rn[0].str.replace("[;#].*$","",regex=True)
            rn[1] = rn[0].str.contains("->")         
            rn[2] = rn[0].str.split("->|=>",expand=True)[0]
            rn[3] = rn[0].str.split("->|=>",expand=True)[1]
            rn[2] = rn[2].str.split("+",expand=False)
            rn[3] = rn[3].str.split("+",expand=False)
        except:
            print("error reading txt file")
        
        # obtaining list of reactions separated as reactive part
        er = list()
        ep = list()
        
        for i in range(len(rn[0])):

            psp = list()
            pst = list()
            ssp = list()
            sst = list()
           
            if rn[2][i]!=['']:
                for j in range(len(rn[2][i])):
                    st=re.search("[0-9]*",rn[2][i][j]).group()
                    if st=='':
                        sst.append(1.0)
                    else:
                        sst.append(float(st))
                    ssp.append(re.search("[^0-9].*",rn[2][i][j]).group())
            else:
                ssp.append([])
                sst.append([])
            
            er.append([ssp,sst])
            
            if rn[3][i]!=['']:
                for j in range(len(rn[3][i])):
                    st=re.search("[0-9]*",rn[3][i][j]).group()
                    if st=='':
                        pst.append(1.0)
                    else:
                        pst.append(float(st))
                    psp.append(re.search("[^0-9].*",rn[3][i][j]).group())
            else:
                psp.append([])
                pst.append([])
            
            ep.append([psp,pst])
            
            if(rn[1][i]):
                er.append([psp,pst])
                ep.append([ssp,sst])
        
        # creation of set of species
        sp=set()
        
        for i in range(len(er)):
           
            if er[i][0][0]:
               r_sp=set(er[i][0])  
            else:
                r_sp=set() 
            
            if ep[i][0][0]:
                p_sp=set(ep[i][0])
            else:
                p_sp=set()
                
            sp=set.union(sp,set.union(r_sp,p_sp))               
        
        # creation of stoichmetric matrix for reactive and product part
        mr=np.zeros([len(er),len(sp)])
        mr=pd.DataFrame(mr,columns=sp)
        mp=copy.copy(mr)
        reac=list()
        prod=list()
        
        for i in range(len(er)):
            r_sp=bt(len(mr.columns))
            r_sp.setall(0)
            
            if er[i][0][0]: 
                for j in range(len(er[i][0])):           
                    ind=mr.columns.get_loc(er[i][0][j])
                    mr.iloc[i][ind]=er[i][1][j]
                    r_sp[ind]=1            
                        
            reac.append(r_sp)
            
            p_sp=bt(len(mr.columns))
            p_sp.setall(0)
            if ep[i][0][0]: 
                for j in range(len(ep[i][0])):
                    ind=mp.columns.get_loc(ep[i][0][j])
                    mp.iloc[i][ind]=ep[i][1][j]
                    p_sp[ind]=1             
 
               
            prod.append(p_sp)              
        # creating class object and retuning it
        out =cls()
        out.sp=mp.columns.values
        out.sp_n=mr.columns.values
        out.mr=mr.T
        out.mp=mp.T
        out.reac=reac
        out.prod=prod
        out.fname=file
        out.txt=True
        out.sbml=False
        
        return out

            
    # Initialization of the network from a smbl file, "file" corresponds 
    # to the path of the smbl file. 
    @classmethod
    def from_sbml(cls,file,modifiers=True,bond_con=True,rand_bc=False):
        try:
            with open(file,'r') as file:
                bs_sbml = bs(file.read(), "xml")
            
            
            er = list()
            ep = list()
            
            
            if bond_con:
                b_i=bs_sbml.findAll(boundaryCondition="true")
            
                if (rand_bc):
                    n_bi=rm.sample(range(len(b_i)),rm.randint(0,len(b_i)))
                else:
                    n_bi=range(len(b_i))
                    
                for j in n_bi:
                    i=b_i[j]
                    er.append([[[]],[[]]])
                    ep.append([[i['id']],[1]])
                
                
            # obtaining list of reactions separated as reactive part 
            # and product part
            bs_r=bs_sbml.select('reaction[id]')            
            # print("bs_r: ",bs_r)
    
            for i in bs_r:

                psp = list()
                pst = list()
                ssp = list()
                sst = list()
                # Calatilic species
                if modifiers:
                    
                    cat = i.find('listOfModifiers')
                    if not (cat is None):
                        cat = cat.select('modifierSpeciesReference')
                        for j in cat:
                            ssp.append(j['species'])
                            sst.append(1)
                            psp.append(j['species'])
                            pst.append(1)
                
                # Reactant species and stoichimetry
                reac = i.find('listOfReactants')
                if not (reac is None):

                    reac = reac.select('speciesReference')
                    for j in reac:

                        ssp.append(j['species'])
                        if j.has_attr('stoichiometry'):
                            sst.append(float(j['stoichiometry'])) 
                        else:
                            sst.append(1.0)
                else:
                        if not ssp:
                            ssp.append([])
                        if not sst:
                            sst.append([])
                
                # Product species and stoichimetry
                prod = i.find('listOfProducts')
                if not (prod is None):
                    prod = prod.select('speciesReference')
                    for j in prod:
                        psp.append(j['species'])
                        if j.has_attr('stoichiometry'):
                            pst.append(float(j['stoichiometry'])) 
                        else:
                            pst.append(1.0)
                else:            
                    if not psp:    
                        psp.append([])
                    if not pst:
                        pst.append([])
                
                er.append([ssp,sst])
                ep.append([psp,pst])
                
                if 'reversible' in i.attrs.keys():
                    if i['reversible'] =='true':
                        ep.append([ssp,sst])
                        er.append([psp,pst])
                        
            # creation of set of species
            sp=set()
            
            for i in range(len(er)):
               
                if er[i][0][0]:
                   r_sp=set(er[i][0])  
                else:
                    r_sp=set() 
                
                if ep[i][0][0]:
                    p_sp=set(ep[i][0])
                else:
                    p_sp=set()
                    
                sp=set.union(sp,set.union(r_sp,p_sp))               
            
            # creation of stoichmetric matrix for reactive and product part
            mr=np.zeros([len(er),len(sp)])
            mr=pd.DataFrame(mr,columns=sp)
            mp=copy.copy(mr)
            reac=list()
            prod=list()
            
            for i in range(len(er)):
                r_sp=bt(len(mr.columns))
                r_sp.setall(0)
                
                if er[i][0][0]: 
                    for j in range(len(er[i][0])):           
                        ind=mr.columns.get_loc(er[i][0][j])
                        mr.iloc[i][ind]=er[i][1][j]
                        r_sp[ind]=1            
                            
                reac.append(r_sp)
                
                p_sp=bt(len(mr.columns))
                p_sp.setall(0)
                if ep[i][0][0]: 
                    for j in range(len(ep[i][0])):
                        ind=mp.columns.get_loc(ep[i][0][j])
                        mp.iloc[i][ind]=ep[i][1][j]
                        p_sp[ind]=1             
     
                   
                prod.append(p_sp)              
            
                 
            # creation of species name vector.
            bs_sp=bs_sbml.select('species[name]')
            sp_n=copy.copy(mp.columns.values)
            for i in bs_sp:
                if i['id'] in mp.columns.values:
                    ind=mp.columns.get_loc(i['id'])
                    sp_n[ind]=i['name']
            
            # creating class object and retuning it
            out=cls()
            out.sp=mp.columns.values
            out.sp_n=sp_n
            out.mr=mr.T
            out.mp=mp.T
            out.prod=prod
            out.reac=reac
            out.fname=file.name
            out.txt=False 
            out.sbml=True 
            
            return out
        except:
            print("error reading smbl file")

    # Function that displays the reactions on the screen. It receives as
    # input a list p of integers corresponding to the reactions to be displayed. 
    # If p is not entered, the complete network is displayed.
    def display(self,p=np.array([])):
        # Creating the random initial rates if it's out of condition
        if len(p)==0:  
            p=self.mp.columns.values
        
        
        for i in p:
            p_text="R_"+str(i)+":   "
            for j in np.where(self.mr.iloc[:,i]!=0)[0]:
                if self.mr.iloc[j,i]==1.0:
                    p_text+=self.mr.index[j]+" "
                elif self.mr.iloc[j,i]==int(self.mr.iloc[j,i]):
                    p_text+=str(int(self.mr.iloc[j,i]))+self.mr.index[j]+" "
                else:
                    p_text+=str(self.mr.iloc[j,i])+self.mr.index[j]+" "
            p_text+="-> "
            for j in np.where(self.mp.iloc[:,i]!=0)[0]:
                if self.mp.iloc[j,i]==1.0:
                    p_text+=self.mp.index[j]+" "
                elif self.mp.iloc[j,i]==int(self.mp.iloc[j,i]):
                    p_text+=str(int(self.mp.iloc[j,i]))+self.mp.index[j]+" "
                else:
                    p_text+=str(self.mp.iloc[j,i])+self.mp.index[j]+" "
                
            print(p_text)
    
    # Funtions that create a mass action dinamics telurrium model (CRNS.model) of 
    # reaction network. It recives as input a vector of the initial concentration 
    # of species i_sp, the reactive constant vector rt and the concentration 
    # thearshold where a species is condidered not present.
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
    
    # def param_model(self, i_sp=np.array([]) ,rt=np.array([])):
    #     # Creating the random initial concetration if it's out of condition
    #     if (i_sp.size==0): 
    #         i_sp=np.random.rand(self.mp.shape[0])
    #     elif(len(i_sp) !=self.mp.shape[0]):
    #         i_sp=np.random.rand(self.mp.shape[0])
        
    #     # Creating the random initial concetration if it's out of condition
    #     if (rt.size==0):  
    #         rt=np.random.rand(self.mp.shape[1])
    #     elif (len(rt) !=self.mp.shape[1]):
    #         rt=np.random.rand(self.mp.shape[1])
         
    #     if self.sbml:
    #         bc_sp=np.where(list(map(lambda x: x in self.model.getBoundarySpeciesIds(), self.mp.index.tolist())))[0]
    #         f_sp=np.where(list(map(lambda x: x in self.model.getFloatingSpeciesIds(), self.mp.index.tolist())))[0]
            
        
    # Function that use a input the existing species and return, the reaction
    # that are be able to be trigered (R_X of X).
    def sp2r(self,sp_set):
        # Checks if input is or not bita array, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set.copy()
        
        nsp=list(set(range(len(sp)))-set(self.bt_ind(sp))) #species no present
        
        rc =[] # creating variable of available reaction 
        for j in range(self.mp.shape[1]):
            if (all(self.mp.iloc[nsp,j]==0) and all(self.mr.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
        
        rbt=bt(self.mp.shape[1])
        rbt.setall(0)
        for i in rc:
            rbt[i]=1
        return(rbt)
    
    
    # Generated closure on a given set, in bitarray or spcies set. "sp_set" 
    # argument can be an np array of species or a bitarray of precent species. 
    # The function  will return an bitarray if bt_type is True, otherwise will
    # return an nparray of species.
    def closure(self,sp_set,bt_type=False):
        
        # Checks if input is or not bita array, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set.copy()
        # creating a vector of reaction that can be trigered
        n_reac = np.array(range(self.mp.shape[1]))
        i=0
        flag=False
        
        # Generates the closure until no reaction can be trigered
        while(len(n_reac)>0):
            if (sp & self.reac[n_reac[i]]) == self.reac[n_reac[i]]:
                
                sp=sp|self.prod[n_reac[i]]
                n_reac=np.delete(n_reac,i)
                flag=True
            else:
                i+=1
            if i>=len(n_reac):
                if(flag):
                    i=0
                    flag=False;
                else:
                    break    
        
        # returns bitarray or species vector
        if bt_type : 
            return sp
        else:
            sp_set=np.array(range(self.mp.shape[0]))
            for i in range(len(sp)):
                if sp[i]==0:
                    sp_set=np.delete(sp_set,np.where(sp_set == i))
            return self.sp[sp_set]
    
    # First iteration of generated closure for a given set, in bitarray or 
    # spcies set. "sp_set" argument can be an np array of species or a bitarray
    # of precent species.The function  will return an bitarray if bt_type is 
    # True, otherwise will return an nparray of species.  
    def clos_one(self,sp_set,bt_type=False):
        
        # Checks if input is or not bita array, if it's no, it make the 
        # transmation
       
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
        # Generates the closure only for the first pass
        for i in range(self.mp.shape[1]):
            if (sp & self.reac[i]) == self.reac[i]:
                sp=sp|self.prod[i]
        
        # returns bitarray or species vector
        if bt_type : 
            return sp
        else:
            sp_set=np.array(range(self.mp.shape[0]))
            for i in range(len(sp)):
                if sp[i]==0:
                    sp_set=np.delete(sp_set,np.where(sp_set == i))
            return self.sp[sp_set]
    
    # Function that confirms if a set is reactive semi-self-mantained, imput is
    # "sp_set" that can be an bitarray or an species array. 
    # returns a true or false depending if the property is achived
    def is_ssm(self,sp_set):
        
        # Checks if input is or not bita array, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
        
        # init of reactant species and product species array
        r = bt(self.mp.shape[0])
        r.setall(0)
        
        p = r.copy()
        p.setall(0)
        
        # generates of reactant species and product species array
        for i in range(self.mp.shape[1]):
             if (sp & self.reac[i]) == self.reac[i]:

                r|=self.reac[i].copy()
                p|=self.prod[i].copy()
        
        # verifies if produces species are in the sp set and if ssm condition 
        # is statisfy

        
        if p.count==0:
            return True;
        elif (sp & p) == p:
            return ((r & p) == r)
        else:
            return False

    # Function that confirms if a set is stoichimetriclly semi-self-mantained, 
    # input is "sp_set" that can be an bitarray or an species array. 
    # returns a true or false depending if the property is achived
    def is_s_ssm(self,sp_set):
        
        # Checks if input is or not bita array, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
            
        else:
            sp=sp_set
              
        # init of reactant species and product species array
        r = bt(self.mp.shape[0])
        r.setall(0)
        
        p = bt(self.mp.shape[0])
        p.setall(0)
        
        # generates of reactant species and product species array
        for i in range(self.mp.shape[1]):
             if (sp & self.reac[i]) == self.reac[i]:
                for j in range(self.mp.shape[0]):
                    d=self.mp.iloc[j][i]-self.mr.iloc[j][i]
                    if d<0:
                        r[j]=1
                    elif d>0:
                        p[j]=1
        
        # verifies if produces species are in the sp set and if ssm condition 
        # is statisfy
        if (sp & p) == p:
            return ((r & p) == r)
        else:
            return False

    def is_sm(self,sp_set):
        
        # Checks if input is or not bitarray, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
            
        # init of reactant species and product species array
        r = bt(self.mp.shape[0])
        r.setall(0)
        
        p = bt(self.mp.shape[0])
        p.setall(0)
        
        # generates of reactant species and product species array
        for i in range(self.mp.shape[1]):
             if (sp & self.reac[i]) == self.reac[i]:
                r|=self.reac[i]
                p|=self.prod[i]
        
        # adding product species to the species vector
        sp=r|p
        
        # init of reaction indexes form trigable reactions
        r_ind=[]
        
        # generates of reaction indexes form trigable reactions
        for i in range(self.mp.shape[1]):
             if (sp & self.reac[i]) == self.reac[i]:
                 r_ind.append(i)
        
        if len(r_ind)==0:
            return True
        
        # init of  indexes form trigable reactions
        sp_ind=[]
        for i in range(len(sp)):
            if sp[i]==1:
                sp_ind.append(i)
        
        # relative Soichiometric matrix to present species and reaction
        S=np.array(self.mr-self.mp)
        S=S[np.ix_(sp_ind,r_ind)]
        S=S.tolist()
        
        c=np.ones(len(r_ind)) #norm of porcess vector for minimization 
        c=c.tolist()
        
        b=np.zeros(len(sp_ind))#imposig that production must be greater or equal to cero
        b=b.tolist()
        
        bounds=[] #creation of bounds vector, solution sould be greater that 1 
        # for each component
        for i in range(len(r_ind)):
            bounds.append((1,None))

        # lineal programing calculation of existance of a solution to be 
        # self-mantainend
        
        res = linprog(c, A_ub=S, b_ub=b, bounds=bounds,method='highs')
        
        return res.success
        
    
    # Subset funcntion for bitarrays
    def is_subset(a, b):
        return (a & b) == a
      
    # Function that returns Bitarray indexes for postition with 1 value 
    def bt_ind(self,bt):
        ind=[]
        for i in range(len(bt)):
            if bt[i]==1:
                ind.append(i)
        
        return ind
    
    
    def gen_basics(self):
        # creating all closed set for clousre of reactant part
        c_reac=[]
          
        for i in self.reac:
            c_reac.append(self.closure(i,True))
        
        # number of equivalence clases
        xeqc=0
        # class equivalence vector
        x_r_p=-np.ones(self.mp.shape[1])
        x_r_p=x_r_p.tolist()
        # list of basic sets
        sp_b=[]
        
        # number of steps
        st=0
        # creating equivalence clases for each reaction that generate the same
        # closure and so each basic can be created
        for i in range(len(x_r_p)):
            st+=1
            # reaction already assigned to equivalnce class
            if x_r_p[i]>=0:
               continue
           
            x_r_p[i]=xeqc
            sp_b.append(c_reac[i])
            
            for j in range(i,self.mp.shape[1]):
                st+=1
                if c_reac[i]==c_reac[j]:
                    x_r_p[j]=xeqc
                    
            xeqc+=1
        
        self.sp_b=sp_b
        self.x_r_p=x_r_p
        
        # assingnig species contained each partition
        sp_p=[]
        
        b_sp=bt(self.mp.shape[0])
        b_sp.setall(0)
        
        for i in range(xeqc):
            sp_p.append(b_sp.copy())
            
        for i in range(self.mp.shape[1]):
            sp_p[x_r_p[i]]|= self.reac[i]|self.prod[i] 
        
        
        self.sp_p=sp_p
        
        
        # assingnig reaction contained each partition
        r_p=[]
        b_r=bt(self.mp.shape[1])
        b_r.setall(0)
        
        for i in range(xeqc):
            for j in np.where(np.array(x_r_p)==i)[0]:
                b_r[j]=1   
            
            r_p.append(b_r.copy())
            b_r.setall(0)
            
            
        self.r_p=r_p
        
        # Creation of other related variables:
        p_b=[] # partitions (equivalence classes) contained in each closure
        r_b=[] # reactions supported by each bassic (equivalence class)
        rsp_b=[] # reactants species contained each bassic (equivalence class)
        psp_b=[] # product species contained each bassic (equivalence class)
        sp_sp_b=[] # stoichiometric positive species contained in the closure of each partition
        sn_sp_b=[] # stoichiometric negative species contained in the closure of each partition
        
        # Assignation of other related variables:
        for i in range(len(self.sp_b)):
            
            p_b.append(self.sp2p(self.sp_b[i]))
            r_b.append(b_r.copy())
            rsp_b.append(b_sp.copy())
            psp_b.append(b_sp.copy())
            sp_sp_b.append(b_sp.copy())
            sn_sp_b.append(b_sp.copy())
            
            for j in self.bt_ind(p_b[i]):
                r_b[i]|=r_p[j]
                for k in self.bt_ind(r_p[j]): 
                    rsp_b[i]|=self.reac[k]
                    psp_b[i]|=self.prod[k]
                    sp_sp_b[i]|=bt(self.mr.iloc[:,k] < self.mp.iloc[:,k])
                    sn_sp_b[i]|=bt(self.mr.iloc[:,k] > self.mp.iloc[:,k])
            
        self.p_b=p_b
        self.r_b=r_b
        self.rsp_b=rsp_b
        self.psp_b=psp_b
        self.sp_sp_b=sp_sp_b
        self.sn_sp_b=sn_sp_b
       
        # Creating connnectig basics and dyncamiclly coected basics sets
        conn=[]
        dyn_conn=[]
        
        # Creating connecting basics and dynamically connected basics sets
        conn=[]
        dyn_conn=[]
        
        b_c=bt(len(self.sp_b))
        b_c.setall(0)
        
        for i in range(len(self.sp_b)):
            
            conn.append(b_c.copy())
            dyn_conn.append(b_c.copy())
            
        for i in range(len(conn)):
            for j in range(len(conn)):
                # connectivity conditions
                if (not j==i) and p_b[i][j]==0 and (rsp_b[i] & psp_b[j]).any():
                    dyn_conn[i][j]=1
                if (not j==i) and p_b[i][j]==0 and (psp_b[i] & rsp_b[j]).any():
                    dyn_conn[j][i]=1
                # Hasse condition (all conected)
                if (not j==i) and p_b[i][j]==0:
                    conn[i][j]=1
                    conn[j][i]=1
        
        self.conn=conn
        self.dyn_conn=dyn_conn
        return(st)
        
    # Function of species that returns the partiton that contains the species
    def sp2p(self, sp):
        p=bt(len(self.sp_b))
        p.setall(0)
        for i in range(len(self.sp_b)):
            if (sp & self.sp_b[i]) == self.sp_b[i]:
                p[i]=1
        
        return p
    
    # Function of partition that returns the contains the species
    def p2sp(self, p):
        sp=bt(len(self.sp))
        sp.setall(0)
        for i in self.bt_ind(p):
            
                sp|=self.sp_b[i]
        
        return sp
    
    #for a given bitset of containend basic set, return the connected basics 
    # sets to it
    def conn_b(self,s):
        
        c=bt(len(self.sp_b))
        c.setall(0)
        
        for i in self.bt_ind(s):
            c|=self.conn[i]
            
        c= c & ~s
        return c
        
    #for a given bitset of containend basic set, return the dynamically 
    # connected basics sets to it
    def dyn_conn_b(self,s):
        c=bt(len(self.sp_b))
        c.setall(0)
        
        for i in self.bt_ind(s):
            c|=self.dyn_conn[i]
            
        c= c & ~s
        return c
    
    #for a given bitset of containend basic set, return the basic sets that 
    # can contribute so the semi-self-maintainace can be reached
    def contrib_b(self,s):
        
        # negative sotichometry species
        
        n=bt(self.mp.shape[0])
        n.setall(0)
        p=n.copy()
        pp=n.copy()
        
        for i in self.bt_ind(s):
            p|=self.sp_sp_b[i]
            n|=self.sn_sp_b[i]
                       
        # P correspond to bitarray whit positive stoichometry
        n= n & ~p
        
        c=s.copy()
        if not n.any():
           c.setall(0)
           return c
       
        # posible sets that can contribute to be semi-self-mantained
        c = ~c
        
        # eliminatig basics form c that no will contrubute
        for i in self.bt_ind(c):
            if not (self.sp_sp_b[i] & n).any():
                c[i]=0
            else:
                pp|=self.sp_sp_b[i]
        
        # pp are posible produces species
        # if they don fullfil consumed species then ther is no contribution
        # to be a semi-self-mantained
        if (n & ~pp).any():
            c.setall(0)
            return c
        else:
            return c
        
    
    # Synergic structure calculation function, requires the gen_basics() 
    # function to be executed beforehand .It returns an directed mulstigraph
    # type networkx (self.syn_str), where each node has as hash the contained basic 
    # bitarray as well attributes as level and set of contained species and if
    # the set si semi-self-maintained. 
    # The edges characterize the closures with the different basic sets, where 
    # the arrival set corresponds to the closure, the key to the 
    # basic with which the closure was performed and the attribute whether the closure 
    # is synergic or not. 
    # The function also retrun the list semi-self-mantained set (self.ssms)
    def gen_syn_str(self):
        if not hasattr(self, 'p_b'):
            print("The basic sets have not been initialized, please run the gen_basics() function.")
            return 
        # Initialization of the synergetic structure as a multigraph object
        G = nx.MultiDiGraph()
        # Initialization of the list of semi-self-maintaine sets
        ssms=[]
        
        #step mesurment
        st=0
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.p_b)):
            st+=1
            G.add_node(fbt(self.p_b[i]),level=self.p_b[i].count(),
                       sp=self.sp[self.bt_ind(self.sp_b[i])],
                       is_ssm=self.is_ssm(self.sp_b[i]))
            if self.is_ssm(self.sp_b[i]):
                ssms.append(self.sp[self.bt_ind(self.sp_b[i])])
        
        # Generation of the multigraph of the synergic structure by set level (contained basics)
        for i in range(len(self.sp_b)):
             
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
             
             # Generating closrues whit connected basics sets for each set in level i
            for j in nodes:
                 for k in self.bt_ind(self.conn_b(bt(j))):
                     st+=1    
                     # Closure result
                     cr_sp=self.closure(self.p2sp(bt(j) | self.p_b[k]),True)
                     cr_p=fbt(self.sp2p(cr_sp))
                     
                     # node is added if si not in structrue
                     if not (cr_p in G):
                         G.add_node(cr_p,level=cr_p.count(),
                                    sp=self.sp[self.bt_ind(cr_sp)],
                                    is_ssm=self.is_ssm(cr_sp))
                         if self.is_ssm(cr_sp):
                             ssms.append(self.sp[self.bt_ind(cr_sp)])
         
                     # Adding edges corresponding to the colsure, and verifing if is a sinergy:
                     if cr_p.count() > (bt(j)|self.p_b[k]).count():
                        G.add_edge(j,cr_p,key=fbt(self.p_b[k]),syn=True)
                     else:
                        G.add_edge(j,cr_p,key=fbt(self.p_b[k]),syn=False)
                
        self.syn_str=G
        self.syn_ssms=ssms
        return(st)
                    

    # Semi-self-maintained structure calculation function, requires the gen_basics() 
    # function to be executed beforehand .It returns an directed multigraph
    # type networkx (self.ssm_str), where each node has as hash the contained basic 
    # bitarray as well attributes as level and set of contained species and if
    # the set si semi-self-maintained. 
    # The edges characterize the closures with the different basic sets, where 
    # the arrival set corresponds to the closure, the key to the 
    # basic with which the closure was performed and the attribute whether the closure 
    # is synergic or not. 
    # The function also retrun the list semi-self-mantained set (self.ssms)
    # This algorithm differs from ge_syn_str() by considering the basic to 
    # be conjugated which can contribute to be semi-self maintained by use of 
    # the contrib_b() function.
    def gen_ssm_str(self):
        if not hasattr(self, 'p_b'):
            print("The basic sets have not been initialized, please run the gen_basics() function.")
            return 
        # Initialization of the synergetic structure as a multigraph object
        G = nx.MultiDiGraph()
        # Initialization of the list of semi-self-maintaine sets
        ssms=[]
        
        # step mesure
        st=0
        
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.p_b)):
            st+=1
            G.add_node(fbt(self.p_b[i]),level=self.p_b[i].count(),
                       sp=self.sp[self.bt_ind(self.sp_b[i])],
                       is_ssm=self.is_ssm(self.sp_b[i]))
            if self.is_ssm(self.sp_b[i]):
                ssms.append(self.sp[self.bt_ind(self.sp_b[i])])
        
        # Generation of the multigraph of the synergic structure by set level (contained basics)
        for i in range(len(self.sp_b)):
             
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
             
             # Generating closrues whit connected basics sets for each set in level i
            for j in nodes:
                # if node is semi-self-maintained (ssm), the it explore other possible 
                # combinations to serch for ssm sets, if not search for basic that can 
                # contibubte to be ssm
                
                if G.nodes[j]["is_ssm"]:
                     conn=self.bt_ind(self.conn_b(bt(j)))
                else:
                     contib=self.contrib_b(bt(j))
                     if not contib.any():
                         continue
                     conn=self.bt_ind(contib)
                 
                for k in conn:
                     st+=1
                     # Closure result
                     cr_sp=self.closure(self.p2sp(bt(j) | self.p_b[k]),True)
                     cr_p=fbt(self.sp2p(cr_sp))
                     
                     # node is added if si not in structrue
                     if not (cr_p in G):
                         G.add_node(cr_p,level=cr_p.count(),
                                    sp=self.sp[self.bt_ind(cr_sp)],
                                    is_ssm=self.is_ssm(cr_sp))
                         if self.is_ssm(cr_sp):
                             ssms.append(self.sp[self.bt_ind(cr_sp)])
         
                     # Adding edges corresponding to the colsure, and verifing if is a sinergy:
                     if cr_p.count() > (bt(j)|self.p_b[k]).count():
                        G.add_edge(j,cr_p,key=k,syn=True)
                     else:
                        G.add_edge(j,cr_p,key=k,syn=False)
                
        self.ssm_str=G
        self.ssm_ssms=ssms
        return(st)
    
    # Semi-self-maintained structure calculation function, requires the gen_basics() 
    # function to be executed beforehand .It returns an directed multigraph
    # type networkx (self.dc_str), where each node has as hash the contained basic 
    # bitarray as well attributes as level and set of contained species and if
    # the set is semi-self-maintained. 
    # The edges characterize the closures with the different basic sets, where 
    # the arrival set corresponds to the closure, the key to the 
    # basic with which the closure was performed and the attribute whether the closure 
    # is synergic or not. 
    # The function also retrun the list semi-self-mantained set (self.ssms)
    # This algorithm differs from ge_syn_str() by considering the basic to 
    # be conjugated which can contribute to be semi-self maintained by use of 
    # the contrib_b() function and use of the function d_connect(), which 
    # only connect to basics if there are reactively connected. The result is a structrure 
    # where the nodes are semi-self-mantianed and only dynamically connected.
    def gen_dyn_str(self):
        if not hasattr(self, 'p_b'):
            print("The basic sets have not been initialized, please run the gen_basics() function.")
            return 
        # Initialization of the synergetic structure as a multigraph object
        G = nx.MultiDiGraph()
        # Initialization of the list of semi-self-maintaine sets
        ssms=[]
        # number of steps
        st=0
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.p_b)):
            st+=1
            G.add_node(fbt(self.p_b[i]),level=self.p_b[i].count(),
                       sp=self.sp[self.bt_ind(self.sp_b[i])],
                       is_ssm=self.is_ssm(self.sp_b[i]))
            if self.is_ssm(self.sp_b[i]):
                ssms.append(self.sp[self.bt_ind(self.sp_b[i])])
        
        # Generation of the multigraph of the synergic structure by set level (contained basics)
        for i in range(len(self.sp_b)):
             
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
             
             # Generating closrues whit connected basics sets for each set in level i
            for j in nodes:
                
                # if node is semi-self-maintained (ssm), the it explore other possible 
                # combinations only whit dynamically conectes basics
                # if not search for basic that can contibubte to be ssm
                
                if G.nodes[j]["is_ssm"]:
                     conn=self.bt_ind(self.dyn_conn_b(bt(j)))
                else:
                     contib=self.contrib_b(bt(j))
                     if not contib.any():
                         continue
                     conn=self.bt_ind(contib)
                 
                for k in conn:
                     st+=1    
                     # Closure result
                     cr_sp=self.closure(self.p2sp(bt(j) | self.p_b[k]),True)
                     cr_p=fbt(self.sp2p(cr_sp))
                     
                     # node is added if si not in structrue
                     if not (cr_p in G):
                         G.add_node(cr_p,level=cr_p.count(),
                                    sp=self.sp[self.bt_ind(cr_sp)],
                                    is_ssm=self.is_ssm(cr_sp))
                         if self.is_ssm(cr_sp):
                             ssms.append(self.sp[self.bt_ind(cr_sp)])
         
                     # Adding edges corresponding to the colsure, and verifing if is a sinergy:
                     if cr_p.count() > (bt(j)|self.p_b[k]).count():
                        G.add_edge(j,cr_p,key=k,syn=True)
                     else:
                        G.add_edge(j,cr_p,key=k,syn=False)
                
        self.dyn_ssm_str=G
        self.dyn_ssms=ssms
        return(st)
    
            
    # Minimal generators generation function, to be started once the basic sets
    # have been generated by the gen_basics() function.  It returns a list
    # (mgen) of bitarray list of species that correspond to the sets so that 
    # by means of the closure they generate the basic set.
    def gen_mgen(self):
        if not hasattr(self, 'p_b'):
            print("The basic sets have not been initialized, please run the gen_basics() function.")
            return 
        
        
        mgen=[]
        
        # Generating a list of the support of each reaction contained in each partition
        for i in range(len(self.r_p)):
            
            v=[]
            for j in self.bt_ind(self.r_p[i]):
                v.append(self.reac[j])
            
            v.sort(key=lambda x: x.count())
        
            # Eliminating supports that contain other supports 
            k=0
            while True: 
                flag=True
                j=k+1
                
                if j >= len(v):
                    break
                while flag:
    
                    if((v[j] & v[k])== v[k]):
                        v.remove(v[j])
                    else:
                        j+=1
                    
                    if not j < len(v):
                        flag=False
                k+=1
                
                if not k < len(v):
                    break
                
            mgen.append(v)
        self.mgen=mgen
    

    # Proto-synergy generation function. it takes as input a minimum generator sp
    # and pi the index of the objective partition. 
    # The function verifies which combinations of partitions can 
    # generate such a generator.        
    def syn_gen(self,sp,pi):

        # list of partition that intersects sp 
        xp=[]     
        
        # union of partition species to se if synergy can be fufill
        ps=bt(len(self.sp))
        ps.setall(0)
        # adding partition that will overlap sp
        for i in range(len(self.sp_b)):
               if (self.sp_p[i] & sp).any():
                   if(i!=pi):
                       xp.append(i)
                       ps|=self.sp_p[i]
        # Verifying if the partitions contain the minimum sp generator 
        # and trigger the proto synergy.
        
        if not ((ps & sp) == sp):
            return
        
        # Bitarray for partition combinations.
        p=bt(len(xp))
        p.setall(0)
        
        # Recursive search of all proto synergies
        for i in range(len(xp)):
            self.r_syn_gen(p,i,sp,xp,pi)
        
            
    # Recursive proto synergy generation function, requires as inputs p the 
    # existing partitions to combine, o the next level to add to the scan, 
    # xp list of indexes of the corresponding partitions, sp minimum 
    # generator to reach and pi the index of the objective partition. 
    # Function recursively explores the possible 
    # combinations to reach a proto synergy. If this is reached, the branch 
    # to be explored will be cut. The proto synergies are stored in a list (syn).
    def r_syn_gen(self,p,o,sp,xp,pi):
        
        # Species bitarray result of the proto sinergy
        u=bt(len(self.sp))
        u.setall(0)
        # Adding partitionn as candiadate to combinate
        p[o]=1;
        
         # Eliminating redundant partitions that are already in account.
        if(p.count()>1):
            ind=self.bt_ind(p).copy()
            for i in ind:
                p[i]=0
                u.setall(0)
                for j in self.bt_ind(p):

                    u|=self.sp_p[xp[j]]
                p[i]=1
            
                if ((self.sp_p[xp[i]] & u & sp) == (self.sp_p[xp[i]] & sp)):
                    p[o]=0
                    return
        
        # Verfing if added partition fulfill triggering the minimal generator
        u|=self.sp_p[xp[o]]
        
        # If it a new synergy, it will be appned and recusrion will stop
        if ((u & sp) == sp):

            c_syn = bt(len(self.p_b))
            c_syn.setall(0)
            
            for i in self.bt_ind(p):

                c_syn[xp[i]]=1
            
            if not c_syn in self.syn:
                self.syn.append(c_syn)
                op=bt(len(self.p_b))
                op.setall(0)
                op[pi]=1
                self.syn_p.append(op)

                
            p[o]=0
            return
        
        # if not recursion continue
        for i in range(o+1,len(xp)):
            self.r_syn_gen(p,i,sp,xp,pi)
     
        p[o]=0
                
            
    # Function that generates all the proto synergies from the minimum 
    # generators. This is achieved through the use of the function gen_syn()
    # and the recursive function r_gen_syn(). The output consists of list syn 
    # which contains all the proto synergies and list syn_p which contains 
    # all the triggered partitions.   
    def all_syn(self):
        if not hasattr(self, 'mgen'):
            print("The minimal genetators have not been initialized, please run the gen_mgen() function.")
            return         
        # List of posible synergies
        
        self.syn=[]
        self.syn_p=[]
        
        # Generation of al proto synergies from all minimum generators
        for i in range(len(self.mgen)):
            for j in self.mgen[i]:
                if j.count()>0:
                    self.syn_gen(j,i)
                else:
                    for k in range(len(self.mgen)):
                        if k!=i:
                            op=bt(len(self.p_b))
                            op.setall(0)
                            op[k]=1
                            self.syn.append(op.copy())
                            op.setall(0)
                            op[i]=1
                            self.syn_p.append(op.copy())
                    
       
    # Verifies which species of a reaction network (not necessarily an organization) 
    # are overproducible. Inputs are species sp_set and process vector pr, 
    # returns a list of overproducible species
    def over_prod(self,sp_set,pr):     
        # Checks if input is or not bita array, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
        
        if not (isinstance(pr,bt)):
            v=bt(self.mp.shape[0])
            v.setall(0)
            
            for i in pr:
                v[i]=1
        else:
            v=pr
        # If it's only one species and self mantained, the reuturns the specie
        # itself
        if sp.count()==1 and self.is_sm(sp):
            return sp
        
        Ns=sp.count() # number of species (rows)
        nsp=list(set(range(len(sp)))-set(self.bt_ind(sp))) #species no present
        
        rc =[] # creating variable of available reaction 
        for j in self.bt_ind(v):
            if (all(self.mp.iloc[nsp,j]==0) and all(self.mr.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
       
        # stoichiometric matrix of rn with prepended destruction reactions
        S=self.mp.iloc[self.bt_ind(sp),rc]-self.mr.iloc[self.bt_ind(sp),rc]
        S=S.to_numpy()
        S=np.column_stack((-np.identity(Ns),S))
        S=S.tolist()
        
        # flow vector constraint: production of every species = 0
        f=np.zeros(Ns) #norm of porcess vector for minimization 
        f=f.tolist()
        
        # cost 0 for every reaction
        cost=np.zeros(Ns+len(rc)) #norm of porcess vector for minimization 
        cost=cost.tolist()
        
        # overproducible status for species (initially False, until proven to be True).
        o=[]
        for i in range(len(sp)):
            o.append(False)
        
        for p in range(Ns):
            if not (o[p]):  # already known overproducible species are skipped 
                
                S[p][p]=1 # creation instead of destruction for p
                f[p]=1 # production of p = 1 instead of 0, if it is possible with creation rate 0, then it is overproducible
                cost[p]=1  # cost for creation of p = 1 instead of 0
            
                # lineal programing calculation of existance of a solution to be 
                # self-mantainend
                
                res = linprog(cost, A_eq=S, b_eq=f,method='highs')
                # The unknown is the process vector v. The equations are S v = f with the inequality constraint v>=0 (rates can be 0).
                # Only the creation reaction for p is penalized in Cost (the rate should be 0 if p is overproducible).
                
                
                if(res.x[p]==0):
                
                    o[p] = True # no need of creation implies p is overproducible
                    for i in range(Ns):
                        if res.x[i]>0:
                            o[i]=True
            
                S[p][p] = -1 
                f[p] = 0
                cost[p] = 0  # original destruction reaction and zero values for next iteration
            
        
        opsp=bt(len(sp))
        opsp.setall(0)
        for i in range(len(o)):
            if o[i]:
                opsp[self.bt_ind(sp)[i]]=1
                
        return(opsp)
    
    # Generates a base of overproduced species of a reaction network 
    # (not necessarily an organization). Inputs are species sp_set and 
    # process vector pr, returns a list of minium overproducible species sets.    
    def op_base(self,sp_set,pr):     
        # Checks if input is or not bita array, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
        
        if not (isinstance(pr,bt)):
            v=bt(self.mp.shape[0])
            v.setall(0)
            
            for i in pr:
                v[i]=1
        else:
            v=pr
            
        # If it's only one species and self mantained, the reuturns the specie
        # itself
        if sp.count()==1 and self.is_sm(sp):
            return {"op_b":[sp]}
        
        Ns=sp.count() # number of species (rows)
        nsp=list(set(range(len(sp)))-set(self.bt_ind(sp))) #species no present
        
        rc =[] # creating variable of available reaction 
        for j in self.bt_ind(v):
            if (all(self.mp.iloc[nsp,j]==0) and all(self.mr.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
       
        # stoichiometric matrix of rn with prepended destruction reactions
        S=self.mp.iloc[self.bt_ind(sp),rc]-self.mr.iloc[self.bt_ind(sp),rc]
        S=S.to_numpy()
        S=np.column_stack((-np.identity(Ns),S))
        S=S.tolist()
        
        
        # flow vector constraint: production of every species = 0
        f=np.zeros(Ns) #norm of porcess vector for minimization 
        f=f.tolist()
        
        # cost 0 for every reaction
        cost=np.zeros(Ns+len(rc)) #norm of porcess vector for minimization 
        cost=cost.tolist()
        
        
        
        #collection for all overproduced bases
        op_b=[]
        
        # coleccion of process vector corresponding to the ovprodued set
        v=[]
        
        for p in range(Ns):
            
            # overproducible status for species (initially False, until proven to be True).
            o=[]
            for i in range(len(sp)):
                o.append(False)
            
            S[p][p]=1 # creation instead of destruction for p
            f[p]=1 # production of p = 1 instead of 0, if it is possible with creation rate 0, then it is overproducible
            cost[p]=1  # cost for creation of p = 1 instead of 0
            
            # lineal programing calculation of existance of a solution to be 
            # self-mantainend
            res = linprog(cost, A_eq=S, b_eq=f,method='highs')
            # The unknown is the process vector v. The equations are S v = f with the inequality constraint v>=0 (rates can be 0).
            # Only the creation reaction for p is penalized in Cost (the rate should be 0 if p is overproducible).

            if(res.x[p]==0):
            
                o[p] = True # no need of creation implies p is overproducible
                for i in range(Ns):
                    if res.x[i]>0:
                        o[i]=True
                
                # adding new overproduced set to the base
                opsp=bt(len(sp))
                opsp.setall(0)
                for i in range(len(o)):
                    if o[i]:
                        opsp[self.bt_ind(sp)[i]]=1
                
                if not (opsp in op_b):
                    op_b.append(opsp)
                    
                    # adding corresponding proces to the base
                    # vi=np.zeros(self.mp.shape[1])
                    # vi[rc]=res.x[Ns:Ns+len(rc)]
                    # v.append(vi)
                    
                
            S[p][p] = -1 
            f[p] = 0
            cost[p] = 0  # original destruction reaction and zero values for next iteration
        
                
        # return({"op_b":op_b,"pr":v})
        return op_b

    # Decomposition function, catalytic species and the respective fragile 
    # cycles of the reaction network. It takes as input the existing 
    # overprduced species opsp, the exiting species sp ande the process vector
    # pr of the present reactions. It returns an array whose components 
    # indicate the  function performed by each species, correlative positions 
    # to sp species  vector. If the value is -1 it corresponds to an 
    # overproducible species, if it is -2 to a catalytic species and if it 
    # is 0 the species is not present. The integer values indicate belonging 
    # to the current fragile cycle.
    def dcom(self,opsp_set,sp_set,pr):
        
        # Checks if input is or not bita array, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
            
        if not (isinstance(opsp_set,bt)):
            opsp=bt(self.mp.shape[0])
            opsp.setall(0)
            
            for i in opsp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    opsp[ind]=1
        else:
            opsp=opsp_set
      
        if not (isinstance(pr,bt)):
            v=bt(self.mp.shape[0])
            v.setall(0)
            
            for i in pr:
                v[i]=1
        else:
            v=pr
        
        # generetinf the decompotion vector
        dcom=np.zeros(len(sp))
        
        # # generating overproducible species
        # opsp = self.over_prod(sp,v)
        dcom[self.bt_ind(opsp)]=-1
        
        # If it's only one species and self mantained, the reuturns the specie
        # itself
        if sp.count()==1 and self.is_sm(sp):
            return dcom
        
        
        # generating matrices to find catalytic and non catalytic species
        sp_ind=self.bt_ind(sp)
        c_m = self.mr.iloc[sp_ind,:].copy()
        nc_m = c_m.copy()
        
        for i in self.bt_ind(v):
            c_m.iloc[:,i] = ((self.mp.iloc[sp_ind,i]!=0) & (self.mr.iloc[sp_ind,i]!=0)) & (self.mp.iloc[sp_ind,i]==self.mr.iloc[sp_ind,i])
            nc_m.iloc[:,i] = (self.mp.iloc[sp_ind,i])!=(self.mr.iloc[sp_ind,i])
        
        
        c_m=c_m.iloc[:,self.bt_ind(v)]
        nc_m=nc_m.iloc[:,self.bt_ind(v)]     
        
        #finding catalytic species
        csp=sp.copy()
        csp.setall(0)
        j=0
        
        
        
        for i in ((c_m.sum(axis=1)!=0) & (nc_m.sum(axis=1)==0)):
            if i:
                csp[sp_ind[j]]=1
            j+=1

        # findinng all other species that aren't catalysers o opverproducble
        fsp=sp.copy() & ~(csp.copy() | opsp.copy())
        
        nsp=list(set(range(len(sp)))-set(self.bt_ind(sp))) #species no present
        
        rc =[] # creating variable of available reaction 
        for j in self.bt_ind(v):
            if (all(self.mp.iloc[nsp,j]==0) and all(self.mr.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
        
        # if there aren't reactive 
        if(len(rc)==0):
            return
        
        # creating fragile cycles
        if fsp.count()>0:            
            # subselecting matrix for fragile cycles species
            m=self.mp.iloc[self.bt_ind(fsp),rc]-self.mr.iloc[self.bt_ind(fsp),rc]
            
            # creating matrix for find conneting species via reactions
            adj=np.zeros((fsp.count(),fsp.count()))
            for i in range(m.shape[1]):
                ka=np.where(m.iloc[:,i]!=0)[0]

                for j in ka:
                    adj[j,ka]=1    
                    
            adj[adj!=0]=1          
            adj = adj + adj.transpose()
            adj[range(adj.shape[0]), range(adj.shape[0])] = 0
            adj=adj>0
            
            fsp_v=np.zeros(fsp.count()) #index of equivalence classes of fragile cycles
            eci=0 
            
            while(any(fsp_v==0)):
                ec = np.where(fsp_v==0)[0][0] # auxilar variable for finding all species of the fragile cycle
                j=0
                while(True):
                    
                    if np.isscalar(ec):
                        ec=np.concatenate((ec,np.where(adj[:,ec])[0]),axis=None)
                    else:
                        ec=np.concatenate((ec,np.where(adj[:,ec[j]])[0]),axis=None)
                        
                    ec=np.unique(ec)
                    j+=1
                    if (j >= len(ec) ):
                        break
                    
                eci+=1 # numbering the equivalance class
                fsp_v[ec]=eci # classifying the corresponding fragile cycle 
           
            # identifying species from the corresponding fragile cycles
            for i in range(fsp.count()):
                dcom[self.bt_ind(fsp)[i]]=fsp_v[i]
                
        # identifying catalyst species
        dcom[self.bt_ind(csp)]=-2
        return dcom
                    
                                       
    # Function that verifiys if a decomposition input is an semi-self-
    # mattained or not.
    def dcom_ssm(self,dcom):
        
        #classification of overproduced species
        opsp=np.where(dcom==-1)[0]

        nsp=np.where(dcom==0)[0] #no present species
        psp=np.where(dcom!=0)[0] #present species
        
        rc =[] # creating variable of available reaction 
        for j in range(self.mp.shape[1]):
            if (all(self.mp.iloc[nsp,j]==0) and all(self.mr.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
        # if there aren't reactive 
        if(len(rc)==0):
            return
        
        # Stochimetric matrix whit available species and reactions
        S=-self.mp.iloc[psp,rc]+self.mr.iloc[psp,rc]
        S=S.values.tolist()
        
        # condition of overproduction for overproduced species
        f=np.zeros(len(psp))
        f[np.isin(psp,opsp)]=-1
        f=f.tolist()
        
        nfcsp=list(set(nsp) | set(opsp)) #non present and opverproduced species, 
        fcrc =[] # creating variable of available reaction for fragile cycles 
        for j in range(self.mp.shape[0]):
            if (all(self.mp.iloc[nfcsp,j]==0) and all(self.mr.iloc[nfcsp,j]==0)):
                fcrc.append(j) # selecting reactions that can be trigger with available species
        # if there aren't reactive 
        # Matrix generation for the equality so fragile cycles species have production 0
        S_eq=(self.mp-self.mr).to_numpy()
        nfcr=list(set(rc)-set(fcrc))
        S_eq[:,nfcr]=0
       
        # cost funtion, all reactions must occur
        c=np.ones(len(rc)).tolist()
        
        bounds=[] #creation of bounds vector, solution sould be greater that 1 
        # for each component
        for i in range(len(rc)):
            bounds.append((1,None))
        
        # lineal programing calculation of existance of a solution to be 
        # self-mantainend
        
        # res = linprog(c, A_ub=S, b_ub=f, A_eq=S_eq, b_eq=f_eq, bounds=bounds)
        res = linprog(c, A_ub=S, b_ub=f, bounds=bounds,method='highs')
        
        # calculation of the production vector
        prod=(-np.array(S).dot(res.x))
        
        # Verifing if overproduced sepcies and fragile cycles species
        if all(prod[f==1]>1) and all(prod[f==0]==0):
            return True
        else:
            return False
                
    # Function that generates a Hasse diagram from a set of species with all 
    # the combinations of the overproduced base, obtaining all the 
    # decompositions. The input corresponds to a set of sp_set species which 
    # must be an organization and the output corresponds to a graph whose 
    # vertexs correspond to the different combinations of overproduced species. 
    # The edges corresponds to the join and meet of the overproduced species. 
    # Each vertex has an attribute that is the decomposition vector. 
    def over_hasse(self,sp_set):
        # Checks if input is or not bitarray, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
      
        # generation of the overproduced base
        op_b=self.op_base(sp_set, self.sp2r(sp_set))
        # initialization of multigraph of opverproduced hasse 
        G = nx.DiGraph()
        all_op=sp.copy()
        all_op.setall(0)
        
        # steps count
        st=0
        
        if len(op_b)==0:
            st+=1
            op=bt(len(sp))
            op.all()
            dcom=self.dcom(op, sp, self.sp2r(sp))
            G.add_node(fbt(op),level=0,
                        dcom=dcom)#,
            return G, st
        
        # The nodes corresponding to the overproduced base.
        for i in op_b:
            st+=1
            dcom=self.dcom(i, sp, self.sp2r(sp))
            G.add_node(fbt(i),level=i.count(),
                        dcom=dcom)#,
                        # is_org=self.dcom_ssm(dcom))
            all_op|=i
        
        

        # print("total levels: ",all_op.count())  
        # Generation of the multigraph of overproduced hasse, by number of 
        # overproduced species
        for i in range(all_op.count()):
            # print("level: ",i)
            # Overproudced sets by level (number of overproduced species)
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
             
            # union whit connected opverproduced bases for each set in level i
            for j in nodes:
                
                # generating connected nodes, nodes that are not containend in j
                conn_nodes=[];
                for l in op_b:
                    if (not (l & bt(j)) == l) and (not fbt(bt(j)|l) in G):
                        conn_nodes.append(l)
                
                for k in conn_nodes:
                    op_new=fbt(bt(j)|k)
                    st+=1
                    # decomposition result
                    dcom=self.dcom(op_new, sp, self.sp2r(sp))
                     
                    # node is added if si not in structrue
                    if not (op_new in G):
                        
                        G.add_node(op_new,level=op_new.count(),
                                  dcom=dcom)#,
                                  # is_org=self.dcom_ssm(dcom))
       
                    # Adding edges corresponding to the hasse diagram
                    if i>0:
                        low_level = [x for x,y in G.nodes(data=True) if y['level']==i-1]
                        for m in low_level:    
                            if bt(m) & bt(op_new) == bt(m):  
                                G.add_edge(m,op_new)
      
        return G, st
        
    # Function that generates a Hasse diagram from a set of species with all 
    # the combinations of the overproduced base, obtaining all the 
    # decompositions. The input corresponds to a set of sp_set species which 
    # must be an organization and the output corresponds to a graph whose 
    # vertexs correspond to the different combinations of overproduced species. 
    # The edges corresponds to the join and meet of the overproduced species. 
    # Each vertex has an attribute that is the decomposition vector. 
    def op_hasse(self,sp_set):
        # Checks if input is or not bitarray, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
      
        # generation of the overproduced base
        op_b=self.op_base(sp_set, self.sp2r(sp_set))
        # initialization of multigraph of opverproduced hasse 
        op_hasse = nx.DiGraph()
        all_op=sp.copy()
        all_op.setall(0)
        
        # The nodes corresponding to the overproduced base.
        for i in op_b:
            dcom=self.dcom(i, sp, self.sp2r(sp))
            op_hasse.add_node(fbt(i),level=i.count(),
                        dcom=dcom)#,
                        # is_org=self.dcom_ssm(dcom))
            all_op|=i
        
        # present opverpruced base elements to combine
        p_op=bt(len(op_b))
        p_op.setall(1)
        
        # Generation of the multigraph of overproduced hasse, by number of 
        # overproduced species
        for i in op_b:
            for j in list(op_hasse.nodes()):
                op_new=i | bt(j)
                # print("exploring: ",op_new)
                if not (op_new in op_hasse):
                    op_hasse.add_node(fbt(op_new),level=op_new.count(),
                                      dcom=self.dcom(op_new, sp, self.sp2r(sp)))#,
                   
        #generating the edges of the graph:
        # for i in range(all_op.count()):
        #     for j in [x for x,y in op_hasse.nodes(data=True) if y['level']==i]:
        #         if i>0:
        #             for m in [x for x,y in op_hasse.nodes(data=True) if y['level']==i-1]:    
        #                 if bt(m) & bt(j) == bt(m):  
        #                     op_hasse.add_edge(m,j)
        print("op_hasse: ", len(op_hasse.nodes()))
        return(op_hasse)

    # Flatten a list using generators comprehensions.
    # Returns a flattened version of list lst.    
    
    def flatten(self,lst: List[Any]) -> Iterable[Any]:


        for sublist in lst:
             if isinstance(sublist, list):
                 for item in sublist:
                     yield item
             else:
                 yield sublist
                 
    # Generates the powerset of an list of elements, returns a list of tuples
    # of each comination posible.
    def powerset(self,iterable):
        s = list(iterable)
        return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
    
    
    # Function that generates a Hasse diagram from a set of species with all 
    # the combinations of the overproduced base, obtaining all the 
    # decompositions. The input corresponds to a set of sp_set species which 
    # must be an organization and the output corresponds to a graph whose 
    # vertexs correspond to the different combinations of overproduced species. 
    # The edges corresponds to the join and meet of the overproduced species. 
    # Each vertex has an attribute that is the decomposition vector. 
    def op_hasse2(self,sp_set):
        # Checks if input is or not bitarray, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
      
        # generation of the overproduced base
        op_b=self.op_base(sp_set, self.sp2r(sp_set))
        # initialization of multigraph of opverproduced hasse 
        op_hasse = nx.DiGraph()
        all_op=sp.copy()
        all_op.setall(0)
        
        # The nodes corresponding to the overproduced base.
        for i in op_b:
            dcom=self.dcom(i, sp, self.sp2r(sp))
            op_hasse.add_node(fbt(i),level=i.count(),
                        dcom=dcom)#,
                        # is_org=self.dcom_ssm(dcom))
            all_op|=i
        
        # generating the powerset and eliminating redudant variables.
        op_b=list(map(self.bt_ind,op_b))
        op_hasse=list(self.powerset(op_b))
        op_hasse=list(map(lambda x: set(list(self.flatten(x))),list(map(list,op_hasse))))
        op_hasse=list(map(sorted,op_hasse))
        op_hasse = [list(x) for x in set(tuple(x) for x in op_hasse)]

        # op_hasse=list(map(list,op_hasse))
        
        # generating all decompositions
        dcom=list(map(lambda x: self.dcom(self.sp[x], sp, self.sp2r(sp)),op_hasse))
        print("op_hasse2: ", len(op_hasse))
        return(op_hasse)
    
    
    
    # Function that generates the polyhedra and polytopes that the fragile and 
    # overproduced cycles of a decomposition. It takes as an increment the 
    # overproduced species opsp_set, the total species sp_set, and the present 
    # reactions pr. It returns as objects the rays and vertices of the 
    # overproduced set op_ph and lists them with the respective rays 
    # and vertices of the fragile circuits fc_ph_list.   
    # def polyth(self,opsp_set,sp_set,pr):
    #     # Checks if input is or not bitarray, if it's no, it make the 
    #     # transmation
    #     if not (isinstance(sp_set,bt)):
    #         sp=bt(self.mp.shape[0])
    #         sp.setall(0)
            
    #         for i in sp_set:
    #             if i in self.mp.index.values:
    #                 ind=self.mp.index.get_loc(i)
    #                 sp[ind]=1
    #     else:
    #         sp=sp_set
        
    #     if not (isinstance(opsp_set,bt)):
    #         opsp=bt(self.mp.shape[0])
    #         opsp.setall(0)
            
    #         for i in opsp_set:
    #             if i in self.mp.index.values:
    #                 ind=self.mp.index.get_loc(i)
    #                 opsp[ind]=1
    #     else:
    #         opsp=opsp_set
            
    #     if not (isinstance(pr,bt)):
    #         v=bt(self.mp.shape[0])
    #         v.setall(0)
            
    #         for i in pr:
    #             v[i]=1
    #     else:
    #         v=pr
        
        
        
    #     Ns=sp.count() # number of species (rows)
    #     nsp=list(set(range(len(sp)))-set(self.bt_ind(sp))) #species no present
        
        
    #     rc =[] # creating variable of available reaction 
    #     for j in self.bt_ind(v):
    #         if (all(self.mp.iloc[nsp,j]==0) and all(self.mr.iloc[nsp,j]==0)):
    #             rc.append(j) # selecting reactions that can be trigger with available species
        
    #     E = np.identity(len(v))
    #     f = np.zeros(len(v))
        
    #     proj =(E,f) 
       
    #     # stoichiometric matrix of rn with prepended destruction reactions
    #     S=self.mp-self.mr
    #     S=S.to_numpy()  
        
    #     S=np.concatenate((S,np.identity(len(v))),axis=0)
        
    #     p=np.zeros(Ns+len(v))
    #     p[self.bt_ind(opsp)]=1
        
    #     ineq=(-S,-p)
        
    #     if set(range(len(v)))-set(rc) != set():
    #         print("d ind: ",list(set(range(len(v)))-set(rc)))
            
    #         d=f.copy()
    #         de=d.copy()
    #         de[np.array(list(set(range(len(v)))-set(rc)))]=1
    #         C=np.diag(de)
            
    #         eq=(C,d)            
    #         op_ph=ph.projection.project_polyhedron(proj, ineq, eq)
    #     else:
            
    #         op_ph=ph.projection.project_polyhedron(proj, ineq)
            
            
    #     dcom=self.dcom(opsp,sp,v)
    #     print("dcom :",dcom)
    #     fc_ph_list=[]
        
    #     for i in np.unique(dcom[dcom>0]):
            
    #         print("i: ",i)
    #         sp_ind=np.where(dcom==i)[0]
    #         # print("sp_ind: ",sp_ind)
            
    #         # print("mr condition: ",self.mr.iloc[sp_ind,:].sum(axis=0)>0)
    #         # print("mp condition: ",self.mp.iloc[sp_ind,:].sum(axis=0)>0)
            
    #         sub_rc=set(np.where(self.mr.iloc[sp_ind,:].sum(axis=0)>0)[0])
    #         sub_rc=set(sub_rc) | set(np.where(self.mp.iloc[sp_ind,:].sum(axis=0))[0])
    #         sub_rc=np.array(list(sub_rc))
   
    #         sub_sp=set(np.where(self.mr.iloc[:,sub_rc].sum(axis=1)>0)[0])
    #         sub_sp=set(sub_sp) | set(np.where(self.mp.iloc[:,sub_rc].sum(axis=1)>0)[0])
    #         # sub_sp=np.array(list(sub_sp))
            
    #         print("sub_sp: ",sub_sp)
            
    #         nsp=list(set(range(len(sp)))-sub_sp) #species no present
    #         sub_rc =[] # creating variable of available reaction 
    #         for j in self.bt_ind(v):
    #             if (all(self.mp.iloc[nsp,j]==0) and all(self.mr.iloc[nsp,j]==0)):
    #                 sub_rc.append(j) # selecting reactions that can be trigger with available species
           
    #         print("sub_rc: ",sub_rc)
            
    #         op_ind=np.array(set(self.bt_ind(opsp)) & set(sp_ind))
    #         cy_ind=np.array(set(dcom[dcom==-2]) & set(sp_ind))
            
    #         E = np.identity(len(v))
    #         f = np.zeros(len(v))
            
    #         proj =(E,f)
        
    #         S=self.mp-self.mr
    #         de=np.zeros(len(v))
            
    #         de[sub_rc]=1
    #         S=np.concatenate((S,np.diag(de)),axis=0)
    #         p=np.zeros(Ns+len(v))
            
    #         # print("S: ",S)
            
    #         ineq=(-S,-p)
            
    #         # condition of free production for overproduced species (either 
    #         # positive or negative)
    #         if set(range(len(v)))-set(sub_rc) != set():
    #             # print("d ind: ",list(set(range(len(v)))-set(sub_rc)))
    #             d=f.copy()
    #             de=d.copy()
    #             de[np.array(list(set(range(len(v)))-set(sub_rc)))]=1
    #             C=np.diag(de)
                    
    #             eq=(C,d)

    #             fc_ph_list.append(ph.projection.project_polyhedron(proj, ineq, eq))
                
    #         else:
    #             fc_ph_list.append(ph.projection.project_polyhedron(proj, ineq))
            
            
        
    #     return({"op_ph":op_ph,"fc_ph_list":fc_ph_list})
            
            
            
            
        
        
    