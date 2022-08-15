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
import matplotlib.pyplot as plt
from pyvis.network import Network
import networkx as nx

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
               # print("a bitarray alone detected, so we are dealing with one set of species will be printed!")
                bit_toList=bit.tolist()
                res="{"
                # Print the species with value 1 in the bit array (1 means present)
                for i in range (len(bit_toList)):
                    if bit_toList[i]==1:
                        res=res+str(self.sp[i])+", "
                #Needed: Delete last ","
                res=res+"}\n"
            elif isinstance(bit,list):
              #  print("a bitarray list detected, so a list of sets of species will be printed")
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
                    res=res
            #Needed: in order to delete last "," we return res[:-3]+"}"
            return(res[:-4]+"}\n")
        except:
           print("input must be \n bitarray representing set of species or \n list of bitarrays representing set of sets of species")
           
    
    # Function that displays the reactions on the screen. It receives as
    # input a list p of integers corresponding to the reactions to be displayed. 
    # If p is not entered, the complete network is displayed.
    def sin_print_r(self,r_set=np.array([])):

        # Checks if input is or not a bitarray, If is, it make the 
        # transformation to an numpy array
        if (r_set.size==0):
            r=self.mp.columns
            r_i=range(len(r))
        else:
            if (isinstance(r_set,bt)):
               r_i=self.bt_ind(r_set)
               r=self.mp.columns[r_i]
            else:
                r_i=[]
                for i in r_set:
                    if i in self.mp.columns:
                        r_i.append(i)
                r=r_set
             
        for i in r_i:
            p_text="R_"+str(i)+":   "
            for j in np.where(self.mr.iloc[:,i]!=0)[0]:
                if self.mr.iloc[j,i]==1.0:
                    p_text+=self.mr.index[j]+" "
                elif self.mr.iloc[j,i]==int(self.mr.iloc[j,i]):
                    p_text+=str(int(self.mr.iloc[j,i]))+self.mr.index[j]+" "
                else:
                    p_text+=str(self.mr.iloc[j,i])+self.mr.index[j]+" "
            p_text+="=> "
            for j in np.where(self.mp.iloc[:,i]!=0)[0]:
                if self.mp.iloc[j,i]==1.0:
                    p_text+=self.mp.index[j]+" "
                elif self.mp.iloc[j,i]==int(self.mp.iloc[j,i]):
                    p_text+=str(int(self.mp.iloc[j,i]))+self.mp.index[j]+" "
                else:
                    p_text+=str(self.mp.iloc[j,i])+self.mp.index[j]+" "
                
            print(p_text)


    # Function that displays the rspecies on the screen. It receives as
    # input a list p of integers or bitarray corresponding to the reactions to be displayed. 
    # If p is not entered, all species is displayed.
    def sin_print_sp(self,sp_set=np.array([])):
        # Checks if input is or not a bitarray, If is, it make the 
        # transformation to an numpy array
        if (sp_set.size==0):
            sp=bt(len(self.sp))
            sp.setall(1)
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set.copy()
        # invoke display function for species
        self.bt_to_sp(sp)
                   
    
    def display_RN(self,r_set=np.array([])):

        # Checks if input is or not a bitarray, If is, it make the 
        # transformation to an numpy array
        if (r_set.size==0):
            r=self.mp.columns
            r_i=range(len(r))
        else:
            if (isinstance(r_set,bt)):
               r_i=self.bt_ind(r_set)
               r=self.mp.columns[r_i]
            else:
                r_i=[]
                for i in r_set:
                    if i in self.mp.columns:
                        r_i.append(i)
                r=r_set
                
        G = nx.MultiDiGraph()
        sp=set()
        for i in r_i:
            sp|=set(self.sp[self.bt_ind(self.reac[i])])|set(self.sp[self.bt_ind(self.prod[i])])
            # size=len(self.sp[self.bt_ind(self.sp_p[i])])*3
            G.add_node("r"+str(i), color = "red", label="r"+str(i), shape="square", size=7)
        
        inf=set(self.is_inflow(np.array(list(sp)),True))
        out=set(self.is_outflow(np.array(list(sp)),True))
        sp-=inf
        sp-=out
        
        for i in sp:
            G.add_node(str(i), color = "blue", label=str(i), size=14, shape="dot")
        for i in inf:
            G.add_node(str(i), color = "green", label=str(i), size=14, shape="dot")
        for i in out:
            G.add_node(str(i), color = "red", label=str(i), size=14, shape="dot")
        for i in r_i:
            
            for j in self.bt_ind(self.reac[i]):
                st_value=self.mr.iloc[j,i]
                label=""
                if st_value!=1:
                    label=str(st_value)
                G.add_edge(str(self.sp[j]), "r"+str(i), color="gray",
                           label=label,title=label)
            
            for j in self.bt_ind(self.prod[i]):
                st_value=self.mp.iloc[j,i]
                label=""
                if st_value!=1:
                    label=str(st_value)
                G.add_edge("r"+str(i), str(self.sp[j]), color="gray",
                           label=label,title=label)
        
        nt = Network('500px', '500px',directed=True)
        nt.from_nx(G)
        nt.toggle_physics(False)
        # nt.show('proto.html')
        return(nt)            
    
    
    #Function that plot the stochimetric matirx, it recives as imput a vector
    # or bit arrar of species (sp_set) and reaction (r_set), and plot the stochimetric 
    # matrix whit colors
    def plot_S(self,sp_set=np.array([]) ,r_set=np.array([])):
        # Checks if input is or not a bitarray, If is, it make the 
        # transformation to an numpy array
        if (sp_set.size==0):
            sp=self.sp
            sp_i=range(len(sp))
        else:
            if (isinstance(sp_set,bt)):
               sp_i=self.bt_ind(sp_set)
               sp=self.sp[sp_i]
            else:
                sp_i=[]
                for i in sp_set:
                    if i in self.mp.index.values:
                        sp_i.append(self.mp.index.get_loc(i))
                sp=sp_set
            
        # Same procedure for reactions
        if (r_set.size==0):
            r=self.mp.columns
            r_i=range(len(r))
        else:
            if (isinstance(r_set,bt)):
               r_i=self.bt_ind(r_set)
               r=self.mp.columns[r_i]
            else:
                r_i=[]
                for i in r_set:
                    if i in self.mp.columns:
                        r_i.append(i)
                r=r_set
            
        # Generating the sotichiometrix sub-matrix        
        S=self.mp.iloc[sp_i,r_i]-self.mr.iloc[sp_i,r_i]
        fig = plt.figure(figsize = (10, 5))
        
        #  Ploting
       
        plt.matshow(S, cmap=plt.cm.viridis)
       
        
        for i in range(len(sp_i)):
            for j in range(len(r_i)):
                c = S.iloc[i,j]               
                plt.text(j, i, str(c), va='center', ha='center')
        
        sp_ticks=np.array(list(map(lambda x: x+ 0.5,range(len(sp_i)))))
        r_ticks=np.array(list(map(lambda x: x+ 0.5,range(len(r_i)))))
       
        plt.yticks(sp_ticks,sp)
        plt.xticks(r_ticks,r)
        plt.show()
            

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
    # Also as input it receives the rpresent eactions r_set as np.array or bitarray format.
    # The function  will return an bitarray if bt_type is True, otherwise will
    # return an numpy.array of species.
    def closure(self,sp_set,r_set=None,bt_type=False):
        
        # Checks if sp_set is or not a bitarray, If not, it make the 
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
            
        # Checks if r_set is or not a bitarray, If not, it make the 
        # transformation to it
        if (r_set is None):
            # creating a vector of reaction that can be triggered
            n_reac = np.array(range(self.mp.shape[1]))
        else:
            if not (isinstance(r_set,bt)):
                r=bt(self.mp.shape[1])
                r.setall(0)
                for i in r_set:
                    r[i]=1
            else:
                r=r_set.copy()

            # creating a vector of reaction that can be triggered
            n_reac = self.bt_ind(r)
        
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
    
    # Function that receives a set of species (sp_set) and returns 
    # the idexes of which are inflow.
    def is_inflow(self,sp_set,set_type=False):
                
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
        
        inf=sp.copy()
        inf.setall(0)
        for i in range(len(self.reac)):
            if (self.reac[i].any()) and (not self.prod[i].any()):
                inf |= self.reac[i]
        
        inf&=sp
        
        if set_type:
            return self.sp[self.bt_ind(inf)]
        
        return self.bt_ind(inf)
        
    # Function that receives a set of species (sp_set) and returns 
    # the idexes of which are inflow.
    def is_outflow(self,sp_set,set_type=False):
                
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
        
        out=sp.copy()
        out.setall(0)
        for i in range(len(self.reac)):
            if (self.prod[i].any()) and (not self.reac[i].any()):
                out |= self.prod[i]
        
        out&=sp
        if set_type:
            return self.sp[self.bt_ind(out)]
        
        return self.bt_ind(out)
    
    # random generation of reaction networks

    # the most simple random network generator, Nr reactions (>1), Ns species (>1)
    # a minimal reaction network is randomly created where each reaction has one reactant and one product and each species
    # is used at least once as reactant and once as product, extra assignments are carried out randomly over reactions
    # there are no inflow or outflow reactions, redundant or null reactions may be generated (rn.merge will filter them out)
    # dist is a log scaled distribution in the [-1,1] range representing locality
    # pr and pp are a log scaled penalization for the repeated use of species as reactants or products
    @classmethod
    def rg_g1(cls,Nr=12,Ns=None,extra=.4, dist=lambda x: x*0+1, pr=0, pp=None):
        
        if Ns is None:
            Ns=Nr
        if pp is None:
            pp=pr
        
        # (0) useful variables and functions
        xr = np.array(list(map(lambda x: x/(Nr-1)*2-1,range(Nr))))  # x coordinates for reactions, range [-1,1]
        xs = np.array(list(map(lambda x: x/(Ns-1)*2-1,range(Ns))))  # x coordinates for species, range [-1,1]
        rnorm = lambda x: x + 2*(x < -1) - 2*(x > 1) # renormalization of coordinates in the range [-1,1]
        norm =lambda x: x / np.sum(x) # Normalization function for probability vector
        usr = np.zeros(Ns)
        usp = usr.copy() # int vectors to count uses of species as reactants and products (initially the counts are 0)
        # (1) creation of the reactants and products Ns X Nr matrices of the reaction network with no species assigned
        sp=list(map(lambda x: "s"+str(x+1),range(Ns)))
        mr=np.zeros([Nr,Ns])
        mr=pd.DataFrame(mr,columns=sp)
        mr=mr.T
        mp=copy.copy(mr)
        # (2) assignment of one random reactant and one random product (different from the reactant) to each reaction
        for i in np.random.choice(range(Nr),Nr,False):
            d=np.array(list(map(dist,list(map(lambda x: rnorm(x-xr[i]),xs)))))
            dr = d - pr*usr; 
            dr = np.exp(dr-np.max(dr))
            sr = np.random.choice(range(Nr),1,p=norm(dr))
            dp = d - pp*usp; 
            dp[sr] = -np.inf 
            dp = np.exp(dp-np.max(dp))
            sp =np.random.choice(range(Ns),1,p=norm(dp)) 
            mr.iloc[sr,i] = 1 
            usr[sr] <- usr[sr]+1
            mp.iloc[sp,i] = 1 
            usp[sp] = usp[sp]+1
        
        # (3) assignment of species not used as reactants to a random reaction (eventually the same reaction)
        i = np.where(usr==0)[0]
        for s in np.random.choice(i,len(i)):
          d=np.array(list(map(dist,list(map(lambda x: rnorm(x-xs[s]),xr)))))
          r = np.random.choice(range(Nr),1,p=norm(np.exp(dr)))
          mr.iloc[s,r] = mr.iloc[s,r] + 1
        
        # (4) assignment of species not used as products to a random reaction (eventually the same reaction)
        i = np.where(usp==0)[0]
        
        for s in np.random.choice(i,len(i)):
          d=np.array(list(map(dist,list(map(lambda x: rnorm(x-xs[s]),xr)))))
          r = np.random.choice(range(Nr),1,p=norm(np.exp(dr)))
          mp.iloc[s,r] = mp.iloc[s,r] + 1
        
        # (5) extra assignment of species at random
        n = np.round(Nr*extra).astype(int) # extra assignments are proportional reactions
        i = np.random.choice(Nr,n,True) # the selected reactions (with repetitions)
        for r in i:
          w = np.random.choice(range(2),2,False)
          d = np.array(list(map(dist,list(map(lambda x: rnorm(x-xr[r]),xs))))) - w[0]*pr*usr - w[1]*pp*usp
          s = np.random.choice(range(Ns),1,p=norm(np.exp(d-np.max(d)))) 
          if w[0]==1: 
              mr.iloc[s,r] = mr.iloc[s,r] + 1
          else: 
              mp.iloc[s,r] = mp.iloc[s,r] + 1
          
          reac=[]
          prod=[]
          for i in range(mp.shape[1]):
              
              r_sp=bt(mr.shape[0])
              r_sp.setall(0)
              p_sp=r_sp.copy()
              
              for j in np.where(mr.iloc[:,i]!=0)[0]:
                  r_sp[j]=1
              for j in np.where(mp.iloc[:,i]!=0)[0]:
                  p_sp[j]=1
              
              reac.append(r_sp)
              prod.append(p_sp)
              
              
          out =cls()
          out.mr=mr
          out.mp=mp
          out.sp=np.array(mr.index)
          out.sp_n=out.sp.copy()
          out.reac=reac
          out.prod=prod
          out.fname=None
          out.txt=False
          out.sbml=False
          
          return out

    
    def rg_extra1(self,p=.1,m=2,Nse=None,extra=None,l="x"):
        
        if Nse is None:
            Nse=np.ceil(self.mr.shape[0]*p).astype(int)
        if extra is None:
            extra=np.round(2*Nse).astype(int)
        
        mr = self.mr
        mp = self.mp
        Ns = mr.shape[0]
        Nr = mr.shape[1]
        # (1) adding extra species
        me=np.zeros((Nr,Nse))
        me=pd.DataFrame(me,columns=list(map(lambda x: l+str(x+1),range(Nse))))
        me=me.T
        mr = pd.concat([mr, me])
        mp = pd.concat([mp, me])
        # (2) using extra species once as reactants and once as products
        for s in range(Nse,(Ns+Nse)):
            mr.iloc[s,np.random.choice(range(Nr),1)] = 1
            mp.iloc[s,np.random.choice(range(Nr),1)] = 1
      
        # (3) extra assignment of species at random
        i = np.random.choice(range(Nr),extra,replace=True) # the selected reactions (with repetitions)
        for r in i:
            s = np.random.choice(range(Nse,(Ns+Nse)),1)
            if np.random.choice(range(2),1)==1: 
                mr.iloc[s,r] = mr.iloc[s,r] + 1
            else:
                mp.iloc[s,r] = mp.iloc[s,r] + 1
        # generation of the new reaction and product bitsets
        reac=[]
        prod=[]
        for i in range(mp.shape[1]):
            
            r_sp=bt(mr.shape[0])
            r_sp.setall(0)
            p_sp=r_sp.copy()
            
            for j in np.where(mr.iloc[:,i]!=0)[0]:
                r_sp[j]=1
            for j in np.where(mp.iloc[:,i]!=0)[0]:
                p_sp[j]=1
            
            reac.append(r_sp)
            prod.append(p_sp)
        
        # Changing the variables of the vectors
        self.mr=mr
        self.mp=mp
        self.sp=np.array(mr.index)
        self.sp_n=self.sp.copy()
        self.reac=reac
        self.prod=prod
    
    # function that adds a percentage of additional (extra) inflow 
    # reactions to an existing network 
    def rg_extra_inflow(self,extra=0.2):
        Ns = self.mr.shape[0]
        
        # selection of species that will considere as inflow species
        i = np.random.choice(range(Ns),np.round(extra*Ns).astype(int),replace=True)
        
        # generating reaction that will be added
        mr=np.zeros((Ns,len(i)))
        mp=mr.copy()
        j=0
        # colname=[]
        for l in i:
            mp[l,j]=1
            # colname.append(Nr+j+1)
            j+=1
            
        mr=pd.DataFrame(mr)
        mr.index=self.mr.index.copy()
        mp=pd.DataFrame(mp)
        mp.index=self.mr.index.copy()
        
        # adding reaction to respective patices
        self.mr = pd.concat([self.mr,mr],axis=1)
        self.mr.columns=range(self.mr.shape[1])
        self.mp = pd.concat([self.mp,mp],axis=1)
        self.mp.columns=range(self.mp.shape[1])
        
        # creating extra bitset variables for the support and products
        reac=[]
        prod=[]
        for i in range(self.mp.shape[1]):
            
            r_sp=bt(self.mr.shape[0])
            r_sp.setall(0)
            p_sp=r_sp.copy()
            
            for j in np.where(self.mr.iloc[:,i]!=0)[0]:
                r_sp[j]=1
            for j in np.where(self.mp.iloc[:,i]!=0)[0]:
                p_sp[j]=1
            
            reac.append(r_sp)
            prod.append(p_sp)       
        
        self.reac=reac
        self.prod=prod
    
    # function that adds a percentage of additional (extra) inflow 
    # reactions to an existing network
    def rg_extra_outflow(self,extra=0.2):
        Ns = self.mr.shape[0]
        
        # selection of species that will considere as outflow species
        i = np.random.choice(range(Ns),np.round(extra*Ns).astype(int),replace=True)
        
        # generating reaction that will be added
        mr=np.zeros((Ns,len(i)))
        mp=mr.copy()
        j=0
        # colname=[]
        for l in i:
            mr[l,j]=1
            # colname.append(Nr+j+1)
            j+=1
            
        mr=pd.DataFrame(mr)
        mr.index=self.mr.index.copy()
        mp=pd.DataFrame(mp)
        mp.index=self.mr.index.copy()
        
        # adding reaction to respective patices
        self.mr = pd.concat([self.mr,mr],axis=1)
        self.mr.columns=range(self.mr.shape[1])
        self.mp = pd.concat([self.mp,mp],axis=1)
        self.mp.columns=range(self.mp.shape[1])
        
        # creating extra bitset variables for the support and products
        reac=[]
        prod=[]
        for i in range(self.mp.shape[1]):
            
            r_sp=bt(self.mr.shape[0])
            r_sp.setall(0)
            p_sp=r_sp.copy()
            
            for j in np.where(self.mr.iloc[:,i]!=0)[0]:
                r_sp[j]=1
            for j in np.where(self.mp.iloc[:,i]!=0)[0]:
                p_sp[j]=1
            
            reac.append(r_sp)
            prod.append(p_sp)       
        
        self.reac=reac
        self.prod=prod

    
    # A wrapper of function rg_g1, but adding percentage of input inflow 
    # and input outflow reactions.
    @classmethod
    def rg_g2(cls,Nr=12,Ns=None,extra=.4, dist=lambda x: x*0+1, pr=0, pp=None, inflow=0.2, outflow=0.2):
        
        # inizialization of 
        out=cls.rg_g1(Nr,Ns,extra, dist, pr, pp)
        out.rg_extra_inflow(inflow)
        out.rg_extra_outflow(outflow)
        out.rn_clean()
        
        return out
            
    # Function that cleans up reaccion redundancies and unsed species.
    def rn_clean(self):
        
        mr=self.mr.copy()
        mp=self.mp.copy()
        
        i = np.where(mr.sum(axis=1)+mr.sum(axis=1)==0)[0]
        if len(i)>0: 
            mr = mr.drop(mr.index[i],axis=0,inplace=False)  # unused species are eliminated
            mp = mp.drop(mp.index[i],axis=0,inplace=False)  # unused species are eliminated
            
        i = np.where(list(map(lambda k: all(mr.iloc[:,k]==mp.iloc[:,k]), range(mr.shape[0]))))[0]
        if len(i)>0:
            mr = mr.drop(mr.index[i],axis=0,inplace=False)  # unused species are eliminated
            mp = mp.drop(mp.index[i],axis=0,inplace=False)  # unused species are eliminated
          
        mr = mr.loc[:,~mr.columns.duplicated()].copy()
        mp = mp.loc[:,~mp.columns.duplicated()].copy()

        # generation of the new reaction and product bitsets
        reac=[]
        prod=[]
        for i in range(mp.shape[1]):
            
            r_sp=bt(mr.shape[0])
            r_sp.setall(0)
            p_sp=r_sp.copy()
            
            for j in np.where(mr.iloc[:,i]!=0)[0]:
                r_sp[j]=1
            for j in np.where(mp.iloc[:,i]!=0)[0]:
                p_sp[j]=1
            
            reac.append(r_sp)
            prod.append(p_sp)
        
        # Changing the variables of the vectors
        self.mr=mr
        self.mp=mp
        self.sp=np.array(mr.index)
        self.sp_n=self.sp.copy()
        self.reac=reac
        self.prod=prod