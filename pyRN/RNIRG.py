#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:52:17 2022

@author: pmaldona

Reaction Network Constructor Class
"""

import numpy as np
import pandas as pd
import re
from bs4 import BeautifulSoup as bs
import copy
from bitarray import bitarray as bt
from scipy.optimize import linprog
import random as rm


# Class definition: It consist to an object generates an stoichiometric matrix 
# for the reactants (RN.mr) as for the products (RN.mp). The default initialization
# is considered form a text file which represents the reaction network as 
# written reaction similar to the antimony format. Also reaction networks 
# can be initialized as smbl
class RNIRG:
    
    # The constructor generates by default the following variables belonging to the object itself:
            # sp -> species array
            # sp_n -> array of the species name (sp_n = sp if not initialized from an sbml file)
            # mr -> stoichiometric matrix of the reactive part
            # mp -> stoichiometric matrix of the productive part
            # reac -> list of the bitarray of species present in the reactive part by reaction
            # prod -> List of the bitarray of species present in the productive part per reaction

    # The following constructors from_txt and from_bml, generate the same 
    # variables from the reading of a file
    def __init__(self):
            self.sp=None
            self.sp_n=None
            self.mr=None
            self.mp=None
            self.prod=None
            self.reac=None

    
    # Reaction Network initialization form text file similar to antimony format, 
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
        
        # creation of stoichiometric matrix for reactive and product part
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
        # creating class object and returning it
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
    # to the path of the smbl file. The variable modifiers=True, includes 
    # these species as catalysts in each reaction, while the condition 
    # bond_cond=True, considers the species with bondarycondition property, 
    # as inflow reactions. Finally, the variable rand_bc=True considers an 
    # increasing number of species with bondarycondition property as inflow. 
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
                
                # Reactant species and stoichiometry
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
                
                # Product species and stoichiometry
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
            
            # creation of stoichiometric matrix for reactive and product part
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
    
    
    #Tomas fork
    #Function that displays the species as a list from a bit array or a list of bitarrays
    #representing a a set of species or a list of sets of species respectively
    #Must be called using print    
    def bt_to_sp(self, bit):
        try:
            # Checks if input is or not a bitarray, If not, it make the 
            # transformation to it
            if isinstance(bit,bt): 
                print("a bitarray alone detected, so we are dealing with one set of species will be printed!")
                bit_toList=bit.tolist()
                res="species set ={"
                # Print the species with value 1 in the bit array (1 means present)
                for i in range (len(bit_toList)):
                    if bit_toList[i]==1:
                        res=res+str(self.sp[i])+", "
                #Needed: Delete last ","
                res=res+"}\n"
            elif isinstance(bit,list):
                print("a bitarray list detected, so a list of sets of species will be printed")
                res=""
                #Apply the process of printing each bitarray in the list
                for i in range(len(bit)):
                    bit_toList=bit[i].tolist()
                    #Here we make reference to the set of species to be printed using 
                    #the index of the list of bitarrays bit
                    res=res+"species set "+str(i)+"={"
                    #Print the species with value 1 in the bit array (1 means present)
                    for j in range (len(bit_toList)):
                        if bit_toList[j]==1:
                            res=res+str(self.sp[j])+", "
                    #Needed: Delete last ","
                    res=res+"}\n"
            return(res)
        except:
           print("input must be \n bitarray representing set of species or \n list of bitarrays representing set of sets of species")
           
    
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

    
    # Function that use a input the existing species and return, the reaction
    # that are be able to be triggered (R_X of X).
    def sp2r(self,sp_set):
        # Checks if input is or not a bitarray, If not, it make the 
        # transformation to it
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
    
    
    # Generated closure of a given set, in bitarray or spices set. "sp_set" 
    # argument can be an numpy.array of species or a bitarray of percent species. 
    # The function  will return an bitarray if bt_type is True, otherwise will
    # return an numpy.array of species.
    def closure(self,sp_set,bt_type=False):
        
        # Checks if input is or not a bitarray, If not, it make the 
        # transformation to it
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set.copy()
        # creating a vector of reaction that can be triggered
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
    # spcies set. "sp_set" argument can be an numpy.array of species or a bitarray
    # of percent species.The function  will return an bitarray if bt_type is 
    # True, otherwise will return an numpy.array of species.  
    def clos_one(self,sp_set,bt_type=False):
        
        # Checks if input is or not a bitarray, If not, it make the 
        # transformation to it
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
            sp_t=sp.copy()
        # Generates the closure only for the first pass
        for i in range(self.mp.shape[1]):
            if (sp_t & self.reac[i]) == self.reac[i]:
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
    
    
    # Function that confirms if a set is reactive semi-self-maintained, input is
    # "sp_set" that can be an bitarray or an species array. 
    # returns a true or false depending if the property is achieved
    def is_ssm(self,sp_set):
        
        # Checks if input is or not a bitarray, If not, it make the 
        # transformation to it
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
        
        # verifies if produces species are in the sp set and if 
        # semi-self-maintained condition is satisfy
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
        
        # Checks if input is or not a bitarray, If not, it make the 
        # transformation to it
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
        
        # Function that confirms if a set is stoichimetriclly semi-self-mantained, 
        # input is "sp_set" that can be an bitarray or an species array. 
        # returns a true or false depending if the property is achived
        if (sp & p) == p:
            return ((r & p) == r)
        else:
            return False
    
    
    # Function that confirms if a set is self-mantained, 
    # input is "sp_set" that can be an bitarray or an species array. 
    # returns a true or false depending if the property is achived
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
