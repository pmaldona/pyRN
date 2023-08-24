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
from bitarray import frozenbitarray as fbt
from scipy.optimize import linprog
import random as rm
import matplotlib.pyplot as plt
from pyvis.network import Network
import networkx as nx
from itertools import combinations
from .genhrn import gen_rn_csp


# Class definition: It consist to an object generates an stoichiometric matrix 
# for the reactants (RN.mr) as for the products (RN.mp). The default initialization
# is considered form a text file which represents the reaction network as 
# written reaction similar to the antimony format. Also reaction networks 
# can be initialized as smbl
class RNIRG:
    
    # The constructor generates by default the following variables belonging to the object itself:
            # SpIdStrArray -> species array
            # SpNameStrArray -> array of the species name (sp_n = sp if not initialized from an sbml file)
            # MrDf -> stoichiometric matrix of the reactive part
            # MpDf -> stoichiometric matrix of the productive part
            # ReacListBt -> list of the bitarray of species present in the reactive part by reaction
            # ProdListBt -> List of the bitarray of species present in the productive part per reaction
            # FilenameStr -> String of the filename, if it correspond
            # IsTextBool -> boolean if is load form a textfile
            # IsSbmlBool -> boolean if its load form a Smbl file

    # The following constructors setFromText and setFromSbml, generate the same 
    # variables reading a file
    def __init__(self):
            self.SpIdStrArray=None
            self.SpNameStrArray=None
            self.MrDf=None
            self.MpDf=None
            self.ProdListBt=None
            self.ReacListBt=None
            self.FilenameStr=None
            self.IsTextBool=False
            self.IsSbmlBool=False
    
    # Reaction Network initialization form text file similar to antimony format, 
    # see example "rn_test.txt". calling the object can be done by RN(file_name). 
    # Example rn = RN.form_txt("rn_text.txt")
    @classmethod
    def setFromText(cls,file):
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
                    st=re.search("^(-?[\d.]+(?:e-?\d+)?)\s*(.*)",rn[2][i][j])
                    try:
                        # st=re.search("[0-9]*",rn[2][i][j]).group()
                        if st is None:
                            sst.append(1.0)
                            ssp.append(re.search("[^0-9].*",rn[2][i][j]).group())
                        else:
                            sst.append(float(st.group(1)))
                            ssp.append(st.group(2))
                        # ssp.append(re.search("[^0-9].*",rn[2][i][j]).group())
                    except:
                        print(st,"0:",st.group(0),"1:",st.group(1),"2:",st.group(2))
            else:
                ssp.append([])
                sst.append([])
            
            er.append([ssp,sst])
            
            if rn[3][i]!=['']:
                for j in range(len(rn[3][i])):
                    st=re.search("^(-?[\d.]+(?:e-?\d+)?)\s*(.*)",rn[3][i][j])
                    try:
                        # st=re.search("[0-9]*",rn[3][i][j]).group()
                        if st is None:
                            pst.append(1.0)
                            psp.append(re.search("[^0-9].*",rn[3][i][j]).group())
                        else:
                            pst.append(float(st.group(1)))
                            psp.append(st.group(2))
                        # psp.append(re.search("[^0-9].*",rn[3][i][j]).group())
                        # psp.append(re.search("(-?[\d.]+(?:e-?\d+)?)\s*(.*)",rn[3][i][j]).group(1))
                    except:
                        print(st,"0:",st.group(0),"1:",st.group(1),"2:",st.group(2))
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
        mr=pd.DataFrame(mr,columns=list(sp))
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
        
        
        #sorting mp and mr dataframes
        sp_i=mp.columns.values.tolist()
        sorted_sp_i=mp.sort_index(axis=1).columns.values.tolist()
        sorted_ind = [sp_i.index(i) for i in sorted_sp_i]
        
        for i in range(len(prod)):
            tmp_prod=bt(len(mr.columns))
            tmp_reac=bt(len(mr.columns))
            for j in range(len(mr.columns)):
                tmp_prod[j]=prod[i][sorted_ind[j]]
                tmp_reac[j]=reac[i][sorted_ind[j]]
            
            prod[i]=tmp_prod
            reac[i]=tmp_reac
        
        mp=mp.sort_index(axis=1)
        mr=mr.sort_index(axis=1)
        
        # creating class object and returning it
        out =cls()
        out.SpIdStrArray=mp.columns.values
        out.SpNameStrArray=mr.columns.values
        out.MrDf=mr.T
        out.MpDf=mp.T
        out.ReacListBt=reac
        out.ProdListBt=prod
        out.FilenameStr=file
        out.IsTextBool=True
        out.IsSbmlBool=False
        
        return out
    
    # Generates a textfile for the reaction network
    def saveToText(self,file,r_set=None):
        
        # Checks if input is or not a bitarray, If is, it make the 
        # transformation to an numpy array
        
        if r_set is not None:
            if (isinstance(r_set,bt)):
               r_i=self.getIndArrayFromBt(r_set)
               r=self.MpDf.columns[r_i]
            else:
                r_i=[]
                for i in r_set:
                    if i in self.MpDf.columns:
                        r_i.append(i)
            
        else:
            r=self.MpDf.columns
            r_i=range(len(r))
        
        self.FilenameStr=file
        text_file = open(self.FilenameStr, "w")
        
        p_text=""
        for i in r_i:
            p_text+="r"+str(i)+":   "
            for j in np.where(self.MrDf.iloc[:,i]!=0)[0]:
                if self.MrDf.iloc[j,i]==1.0:
                    p_text+=self.MrDf.index[j]+" "
                elif self.MrDf.iloc[j,i]==int(self.MrDf.iloc[j,i]):
                    p_text+=str(int(self.MrDf.iloc[j,i]))+self.MrDf.index[j]+" "
                else:
                    p_text+=str(self.MrDf.iloc[j,i])+self.MrDf.index[j]+" "
                p_text+="+ "
            if len(np.where(self.MrDf.iloc[:,i]!=0)[0])>0:
                p_text=p_text[:-2]
            p_text+="=> "
            
            for j in np.where(self.MpDf.iloc[:,i]!=0)[0]:
                if self.MpDf.iloc[j,i]==1.0:
                    p_text+=self.MpDf.index[j]+" "
                elif self.MpDf.iloc[j,i]==int(self.MpDf.iloc[j,i]):
                    p_text+=str(int(self.MpDf.iloc[j,i]))+self.MpDf.index[j]+" "
                else:
                    p_text+=str(self.MpDf.iloc[j,i])+self.MpDf.index[j]+" "
                p_text+="+ "
            if len(np.where(self.MpDf.iloc[:,i]!=0)[0])>0:
                p_text=p_text[:-2]
            
            p_text+="\n"
        text_file.write(p_text)
        text_file.close()
            
    # Initialization of the network from a smbl file, "file" corresponds 
    # to the path of the smbl file. The variable modifiers=True, includes 
    # these species as catalysts in each reaction, while the condition 
    # bond_cond=True, considers the species with bondarycondition property, 
    # as inflow reactions. Finally, the variable rand_bc=True considers an 
    # increasing number of species with bondarycondition property as inflow. 
    @classmethod
    def setFromSbml(cls,file,modifiers=True,bond_con=True,rand_bc=False):
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
            mr=pd.DataFrame(mr,columns=list(sp))
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
            
            #sorting mp and mr dataframes
            sp_i=mp.columns.values.tolist()
            sorted_sp_i=mp.sort_index(axis=1).columns.values.tolist()
            sorted_ind = [sp_i.index(i) for i in sorted_sp_i]
            
            for i in range(len(prod)):
                tmp_prod=bt(len(mr.columns))
                tmp_reac=bt(len(mr.columns))
                for j in range(len(mr.columns)):
                    tmp_prod[j]=prod[i][sorted_ind[j]]
                    tmp_reac[j]=reac[i][sorted_ind[j]]
                
                prod[i]=tmp_prod
                reac[i]=tmp_reac
            
            mp=mp.sort_index(axis=1)
            mr=mr.sort_index(axis=1)
                 
            # creation of species name vector.
            bs_sp=bs_sbml.select('species[name]')
            sp_n=copy.copy(mp.columns.values)
            for i in bs_sp:
                if i['id'] in mp.columns.values:
                    ind=mp.columns.get_loc(i['id'])
                    sp_n[ind]=i['name']
            
   
            
            # creating class object and retuning it
            out=cls()
            out.SpIdStrArray=mp.columns.values
            out.SpNameStrArray=sp_n
            out.MrDf=mr.T
            out.MpDf=mp.T
            out.ProdListBt=prod
            out.ReacListBt=reac
            out.FilenameStr=file.name
            out.IsTextBool=False 
            out.IsSbmlBool=True 
            
            return out
        except:
            print("error reading smbl file")
    
    def setSpConnMat(self,transitive=True):
        '''
        

        Returns
        -------
        Generates a connectivity graph of species via reactions.

        '''
        
        #creation species coneted graph
        self.SpRConnNx=nx.Graph()
        # adding nodes to the graph
        for i in self.SpIdStrArray:
            self.SpRConnNx.add_node(i)
        
        #generating directly connected species in a reaction
        for i in range(len(self.ReacListBt)):
            conn_sp=self.ReacListBt[i] | self.ProdListBt[i]     
            sp_indx=conn_sp.search(1)
            
            for j in range(len(sp_indx)):
                b_sp=self.SpIdStrArray[sp_indx[j]]
                for k in sp_indx[j:]:
                    f_sp=self.SpIdStrArray[k]
                    if not self.SpRConnNx.has_edge(b_sp,f_sp):
                        self.SpRConnNx.add_edge(b_sp,f_sp,w=1)
                    else:                        
                        self.SpRConnNx.edges[b_sp, f_sp]["w"]+=1
        self.SpConnNx=nx.transitive_closure(self.SpRConnNx,reflexive=True)
            
    def setSpCConnMat(self):
        '''
        

        Returns
        -------
        Generates a causal connectivity graph of species via reactions.

        '''
        
        #creation species coneted graph
        self.SpCConnNx=nx.DiGraph()
        # adding nodes to the graph
        for i in self.SpIdStrArray:
            self.SpCConnNx.add_node(i)
        
        #generating directly connected species in a reaction
        for i in range(len(self.ReacListBt)):
            f_sp=self.ReacListBt[i]      
            t_sp=self.ProdListBt[i]
            
            
            for j in f_sp.search(1):
                for k in t_sp.search(1):
                    f_sp=self.SpIdStrArray[k]
                    if not self.SpCConnNx.has_edge(self.SpIdStrArray[j],self.SpIdStrArray[k]):
                        self.SpCConnNx.add_edge(self.SpIdStrArray[j],self.SpIdStrArray[k],w=1)
                    else:                        
                        self.SpCConnNx.edges[self.SpIdStrArray[j], self.SpIdStrArray[k]]["w"]+=1
        
    
    def setSpConnFrac(self):
        '''
            
    
        Returns
        -------
        Generates fraction of posible connectivity for species.
    
        '''
        
        # total edges
        total_edges=self.SpConnNx.number_of_nodes() * (self.SpConnNx.number_of_nodes() + 1) / 2
        
        # Conneted fration
        self.SpConnFrac = len(self.SpConnNx.edges) / total_edges

    
    def getConnSp(self,sp_set,bitout=False):
        '''
        

        Parameters
        ----------
        sp_set : bitarray ot SpId array
            Candiadate species we want to explore.
        bitout: bool,
            If True, then it returns a bitarra of the present sepcies, if False, 
            it return a numpty arrat of the species id.

        Returns
        -------
        Connected species to sp_set.

        '''
        #Check type of file and set it to species array
        if sp_set is None:
            sp_set=self.SpIdStrArray
            
        if (isinstance(sp_set,bt)):
            bit_list=sp_set.tolist()
            sp=[]
            for i in range(len(bit_list)):
                if bit_list[i]==1:
                    sp.append(self.MrDf.index[i])
            # sp=np.array(sp)
        
        # print(sp)
        # Add variables in consideration of the species connection graph 
        out=set()
        for i in sp:
            # print(i)
            # print(set(dict(self.SpConnNx['s1']).keys()))
            out|=set(dict(self.SpConnNx[i]).keys())
        
        out-=set(sp)
        out=np.array(list(out))
        
        # Transform the output fo bitarray if asked
        if bitout:
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in out:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
            out=sp
        
        return(out)
  
    
    
    
    #Tomas fork
    #Function that displays the species as a list from a bit array or a list of bitarrays
    #representing a a set of species or a list of sets of species respectively
    #Must be called using print    
    def printSpIdFromBt(self, bit):
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
                        res=res+str(self.SpIdStrArray[i])+", "
                #Needed: Delete last ","
                res=res[:-1]+"}"
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
                            res=res+str(self.SpIdStrArray[j])+", "
                    res=res[:-2]
                    res=res+"} "
            #Needed: in order to delete last "," we return res[:-3]+"}"
            return(res)
        except:
           print("input must be \n bitarray representing set of species or \n list of bitarrays representing set of sets of species")
           
    
    # Function that displays the reactions on the screen. It receives as
    # input a list p of integers corresponding to the reactions to be displayed. 
    # If p is not entered, the complete network is displayed.
    def printRp(self,r_set=None,string_out=False):

        # Checks if input is or not a bitarray, If is, it make the 
        # transformation to an numpy array

        if (isinstance(r_set,bt)):
           r_i=self.getIndArrayFromBt(r_set)
           # r=self.MpDf.columns[r_i]
        else:
            if r_set is None:
                # r=self.MpDf.columns
                r_i=range(self.MpDf.shape[1])
            else:
                r_i=[]
                for i in r_set:
                    if i in self.MpDf.columns:
                        r_i.append(i)
                # r=r_set
            
        str_out=str()
        for i in r_i:
            p_text="r"+str(i)+":   "
            for j in np.where(self.MrDf.iloc[:,i]!=0)[0]:
                if self.MrDf.iloc[j,i]==1.0:
                    p_text+=self.MrDf.index[j]+" "
                elif self.MrDf.iloc[j,i]==int(self.MrDf.iloc[j,i]):
                    p_text+=str(int(self.MrDf.iloc[j,i]))+self.MrDf.index[j]+" "
                else:
                    p_text+=str(self.MrDf.iloc[j,i])+self.MrDf.index[j]+" "
                p_text+="+ "
            if len(np.where(self.MrDf.iloc[:,i]!=0)[0])>0:
                p_text=p_text[:-2]
            p_text+="=> "
            
            for j in np.where(self.MpDf.iloc[:,i]!=0)[0]:
                if self.MpDf.iloc[j,i]==1.0:
                    p_text+=self.MpDf.index[j]+" "
                elif self.MpDf.iloc[j,i]==int(self.MpDf.iloc[j,i]):
                    p_text+=str(int(self.MpDf.iloc[j,i]))+self.MpDf.index[j]+" "
                else:
                    p_text+=str(self.MpDf.iloc[j,i])+self.MpDf.index[j]+" "
                p_text+="+ "
            if len(np.where(self.MpDf.iloc[:,i]!=0)[0])>0:
                p_text=p_text[:-2]    
            if not string_out:
                print(p_text)
            else:
                str_out+=p_text+"\n"
            
        if string_out:
            return str_out
        
    # Function that displays the reactions on the screen. It receives as
    # input a list p of integers corresponding to the reactions to be displayed. 
    # If p is not entered, the complete network is displayed.
    def printRpFromProcess(self,r_set=np.array([]),pr=np.array([])):
        
        if (isinstance(r_set,bt)):
           r_i=self.getIndArrayFromBt(r_set)
           # r=self.MpDf.columns[r_i]
        else:
            if (r_set.size==0):
                # r=self.MpDf.columns
                r_i=range(self.MpDf.shape[1])
            else:
                r_i=[]
                for i in r_set:
                    if i in self.MpDf.columns:
                        r_i.append(i)
                # r=r_set
                 
        # Checks if input is or not a 
        for i in r_i:
            p_text="r"+str(i)+":   "
            for j in np.where(self.MrDf.iloc[:,i]!=0)[0]:
                if self.MrDf.iloc[j,i]==1.0:
                    p_text+=self.MrDf.index[j]+" "
                elif self.MrDf.iloc[j,i]==int(self.MrDf.iloc[j,i]):
                    p_text+=str(int(self.MrDf.iloc[j,i]))+self.MrDf.index[j]+" "
                else:
                    p_text+=str(self.MrDf.iloc[j,i])+self.MrDf.index[j]+" "
                p_text+="+ "
            if len(np.where(self.MrDf.iloc[:,i]!=0)[0])>0:
                p_text=p_text[:-2]
            p_text+="("+str(pr[i])+") => "
            
            for j in np.where(self.MpDf.iloc[:,i]!=0)[0]:
                if self.MpDf.iloc[j,i]==1.0:
                    p_text+=self.MpDf.index[j]+" "
                elif self.MpDf.iloc[j,i]==int(self.MpDf.iloc[j,i]):
                    p_text+=str(int(self.MpDf.iloc[j,i]))+self.MpDf.index[j]+" "
                else:
                    p_text+=str(self.MpDf.iloc[j,i])+self.MpDf.index[j]+" "
                p_text+="+ "
            if len(np.where(self.MpDf.iloc[:,i]!=0)[0])>0:
                p_text=p_text[:-2]    
            print(p_text)

    
    # Function that displays the rspecies on the screen. It receives as
    # input a list p of integers or bitarray corresponding to the reactions to be displayed. 
    # If p is not entered, all species is displayed.
    def printSp(self,sp_set=None):
        # Checks if input is or not a bitarray, If is, it make the 
        # transformation to an numpy array
        if sp_set is None:
            sp_set=bt(len(self.SpIdStrArray))
            sp_set.setall(1)
            
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set.copy()
        # invoke display function for species
        print("Species: ",self.printSpIdFromBt(sp))
        
                   
    
    def getRnDisplayPv(self,r_set=None,x_size='500px',y_size='500px',notebook=False,cdn_resources='local'):

        # Checks if input is or not a bitarray, If is, it make the 
        # transformation to an numpy array
        if r_set is None:
            r=self.MpDf.columns
            r_i=range(len(r))
        else:
            if (isinstance(r_set,bt)):
               r_i=self.getIndArrayFromBt(r_set)
               r=self.MpDf.columns[r_i]
            else:
                r_i=[]
                for i in r_set:
                    if i in self.MpDf.columns:
                        r_i.append(i)
                r=r_set
                
        G = nx.MultiDiGraph()
        sp=set()
        for i in r_i:
            sp|=set(self.SpIdStrArray[self.getIndArrayFromBt(self.ReacListBt[i])])|set(self.SpIdStrArray[self.getIndArrayFromBt(self.ProdListBt[i])])
            # size=len(self.SpIdStrArray[self.getIndArrayFromBt(self.SpIdStrArray_p[i])])*3
            G.add_node("r"+str(i), color = "yellow", label="r"+str(i), shape="square", size=7)
        
        
        not_r_i=set(range(self.MpDf.shape[1]))-set(r_i)
        for i in not_r_i:
            G.add_node("r"+str(i), color = "#E0E0E0", label="r"+str(i), shape="square", size=7)
        
        
        nsp=set(self.SpIdStrArray)-sp
        
        
        
        inf=set(self.getInflowFromSp(np.array(list(sp)),"id"))
        out=set(self.getOutflowFromSp(np.array(list(sp)),"id"))
        sp-=inf
        sp-=out
        
        for i in sp:
            G.add_node(str(i), color = "blue", label=str(i), size=14, shape="dot")
        for i in nsp:
            G.add_node(str(i), color = "#E0E0E0", label=str(i), size=14, shape="dot")
        for i in inf:
            G.add_node(str(i), color = "green", label=str(i), size=14, shape="dot")
        for i in out:
            G.add_node(str(i), color = "red", label=str(i), size=14, shape="dot")
        for i in r_i:
            
            for j in self.getIndArrayFromBt(self.ReacListBt[i]):
                st_value=self.MrDf.iloc[j,i]
                label=""
                if st_value!=1:
                    label=str(st_value)
                G.add_edge(str(self.SpIdStrArray[j]), "r"+str(i), color="gray",
                           label=label,title=label)
            
            for j in self.getIndArrayFromBt(self.ProdListBt[i]):
                st_value=self.MpDf.iloc[j,i]
                label=""
                if st_value!=1:
                    label=str(st_value)
                G.add_edge("r"+str(i), str(self.SpIdStrArray[j]), color="gray",
                           label=label,title=label)
        
        for i in not_r_i:
            
            for j in self.getIndArrayFromBt(self.ReacListBt[i]):
                st_value=self.MrDf.iloc[j,i]
                label=""
                if st_value!=1:
                    label=str(st_value)
                G.add_edge(str(self.SpIdStrArray[j]), "r"+str(i), color="#E0E0E0",
                           label=label,title=label)
            
            for j in self.getIndArrayFromBt(self.ProdListBt[i]):
                st_value=self.MpDf.iloc[j,i]
                label=""
                if st_value!=1:
                    label=str(st_value)
                G.add_edge("r"+str(i), str(self.SpIdStrArray[j]), color="#E0E0E0",
                           label=label,title=label)
            
        
        nt = Network(x_size, y_size ,directed=True,notebook=notebook,cdn_resources=cdn_resources)
        nt.from_nx(G)
        nt.toggle_physics(False)
        # nt.show('proto.html')
        return(nt)            
    
    
    #Function that plot the stochimetric matirx, it recives as imput a vector
    # or bit arrar of species (sp_set) and reaction (r_set), and plot the stochimetric 
    # matrix With colors
    def plotS(self,sp_set=np.array([]) ,r_set=np.array([]), return_figure=False):
        # Checks if input is or not a bitarray, If is, it make the 
        # transformation to an numpy array
        if (sp_set.size==0):
            sp=self.SpIdStrArray
            sp_i=range(len(sp))
        else:
            if (isinstance(sp_set,bt)):
               sp_i=self.getIndArrayFromBt(sp_set)
               sp=self.SpIdStrArray[sp_i]
            else:
                sp_i=[]
                for i in sp_set:
                    if i in self.MpDf.index.values:
                        sp_i.append(self.MpDf.index.get_loc(i))
                sp=sp_set
            
        # Same procedure for reactions
        if (r_set.size==0):
            r=self.MpDf.columns
            r_i=range(len(r))
        else:
            if (isinstance(r_set,bt)):
               r_i=self.getIndArrayFromBt(r_set)
               r=self.MpDf.columns[r_i]
            else:
                r_i=[]
                for i in r_set:
                    if i in self.MpDf.columns:
                        r_i.append(i)
                r=r_set
            
        # Generating the sotichiometrix sub-matrix        
        S=self.MpDf.iloc[sp_i,r_i]-self.MrDf.iloc[sp_i,r_i]
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
        
        # Needed for plotting matplotlab figures in the UI
        if return_figure == True:
            return plt
        else:
            plt.show()
            

    # Function that use a input the existing species and return, the reaction
    # that are be able to be triggered (R_X of X).
    def getRpFromSp(self,sp_set):
        # Checks if input is or not a bitarray, If not, it make the 
        # transformation to it
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set.copy()
        
        nsp=list(set(range(len(sp)))-set(self.getIndArrayFromBt(sp))) #not present species
        
        rc =[] # creating variable of available reaction 
        for j in range(self.MpDf.shape[1]):
            if (all(self.MpDf.iloc[nsp,j]==0) and all(self.MrDf.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
        
        rbt=bt(self.MpDf.shape[1])
        rbt.setall(0)
        for i in rc:
            rbt[i]=1
        return(rbt)
    
    def getTriggerableRpBtFromSp(self,sp_set):
        '''
        

        Parameters
        ----------
        sp_set : bitset or numpy array
            Bitset of species of nupy array of present species.

        Returns
        -------
        A bitarray of triggerable reaction by sp_set species set.

        '''
        
        # Checks if input is or not a bitarray, If not, it make the 
        # transformation to it
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set.copy()
        
        nsp=list(set(range(len(sp)))-set(self.getIndArrayFromBt(sp))) #not present species
        
        rc =[] # creating variable of available reaction 
        for j in range(self.MpDf.shape[1]):
            if (all(self.MrDf.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
        
        rbt=bt(self.MpDf.shape[1])
        rbt.setall(0)
        for i in rc:
            rbt[i]=1
        return(rbt)
    
    # Generated closure of a given set, in bitarray or spices set. "sp_set" 
    # argument can be an numpy.array of species or a bitarray of percent species. 
    # Also as input it receives the rpresent eactions r_set as np.array or bitarray format.
    # The function  will return an bitarray if bt_type is True, otherwise will
    # return an numpy.array of species.
    def getClosureFromSp(self,sp_set,r_set=None,bt_type=False):
        
        # Checks if sp_set is or not a bitarray, If not, it make the 
        # transformation to it
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set.copy()
            
        # Checks if r_set is or not a bitarray, If not, it make the 
        # transformation to it
        if (r_set is None):
            # creating a vector of reaction that can be triggered
            n_reac = np.array(range(self.MpDf.shape[1]))
        else:
            if not (isinstance(r_set,bt)):
                r=bt(self.MpDf.shape[1])
                r.setall(0)
                for i in r_set:
                    r[i]=1
            else:
                r=r_set.copy()

            # creating a vector of reaction that can be triggered
            n_reac = self.getIndArrayFromBt(r)
        
        i=0
        flag=False
        # Generates the closure until no reaction can be trigered
        while(len(n_reac)>0):
            if (sp & self.ReacListBt[n_reac[i]]) == self.ReacListBt[n_reac[i]]:
                
                sp=sp|self.ProdListBt[n_reac[i]]
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
            sp_set=np.array(range(self.MpDf.shape[0]))
            for i in range(len(sp)):
                if sp[i]==0:
                    sp_set=np.delete(sp_set,np.where(sp_set == i))
            return self.SpIdStrArray[sp_set]
    
    
    # First iteration of generated closure for a given set, in bitarray or 
    # spcies set. "sp_set" argument can be an numpy.array of species or a bitarray
    # of percent species.The function  will return an bitarray if bt_type is 
    # True, otherwise will return an numpy.array of species.  
    def getClosOneFromSp(self,sp_set,bt_type=False):
        
        # Checks if input is or not a bitarray, If not, it make the 
        # transformation to it
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
            sp_t=sp.copy()
        # Generates the closure only for the first pass
        for i in range(self.MpDf.shape[1]):
            if (sp_t & self.ReacListBt[i]) == self.ReacListBt[i]:
                sp=sp|self.ProdListBt[i]
        
        # returns bitarray or species vector
        if bt_type : 
            return sp
        else:
            sp_set=np.array(range(self.MpDf.shape[0]))
            for i in range(len(sp)):
                if sp[i]==0:
                    sp_set=np.delete(sp_set,np.where(sp_set == i))
            return self.SpIdStrArray[sp_set]
    
    def getNonReacSets(self,sp_set):
        '''
        

        Parameters
        ----------
        sp_set : bitarray
            set of species.

        Returns
        -------
        A list of all bitarray sets of species, whose additive part is 
        non-reactive with respect to sp_set.

        '''
        
        non_org_reac=~self.getTriggerableRpBtFromSp(sp_set)
        non_org=~sp_set
        
        cand_sets = [com for sub in range(1,non_org.count()) for com in combinations(self.getIndArrayFromBt(non_org), sub)]
        cand_sets = list(map(lambda x: fbt(self.getBtFromIndArray(x,len(sp_set))),cand_sets))
        
        non_set=list(map(lambda x: fbt(self.ReacListBt[x]),non_org_reac.search(1)))
        
        out=cand_sets.copy()
        for i in cand_sets:
            if any(list(map(lambda x: (x & i) == x ,non_set))):
                out.remove(i)

        return(list(map(lambda x: bt(x|sp_set),out)))
    
    
    # Function that confirms if a set is reactive semi-self-maintained, input is
    # "sp_set" that can be an bitarray or an species array. 
    # returns a true or false depending if the property is achieved
    def isSsmFromSp(self,sp_set):
        
        # Checks if input is or not a bitarray, If not, it make the 
        # transformation to it
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
        
        # init of reactant species and product species array
        r = bt(self.MpDf.shape[0])
        r.setall(0)
        
        p = r.copy()
        p.setall(0)
        
        # generates of reactant species and product species array
        for i in range(self.MpDf.shape[1]):
             if (sp & self.ReacListBt[i]) == self.ReacListBt[i]:
                 
                r|=self.ReacListBt[i].copy()
                p|=self.ProdListBt[i].copy()
        
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
    def isStoiSsmFromSp(self,sp_set):
        
        # Checks if input is or not a bitarray, If not, it make the 
        # transformation to it
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
            
        else:
            sp=sp_set
              
        # init of reactant species and product species array
        r = bt(self.MpDf.shape[0])
        r.setall(0)
        
        p = bt(self.MpDf.shape[0])
        p.setall(0)
        
        # generates of reactant species and product species array
        for i in range(self.MpDf.shape[1]):
             if (sp & self.ReacListBt[i]) == self.ReacListBt[i]:
                for j in range(self.MpDf.shape[0]):
                    d=self.MpDf.iloc[j][i]-self.MrDf.iloc[j][i]
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
    def isSmFromSp(self,sp_set,ret_vect=False):
        
        # Checks if input is or not bitarray, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
            
        # init of reactant species and product species array
        r = bt(self.MpDf.shape[0])
        r.setall(0)
        
        p = bt(self.MpDf.shape[0])
        p.setall(0)
        
        # generates of reactant species and product species array
        for i in range(self.MpDf.shape[1]):
             if (sp & self.ReacListBt[i]) == self.ReacListBt[i]:
                r|=self.ReacListBt[i]
                p|=self.ProdListBt[i]
        
        # adding product species to the species vector
        sp=r|p
        
        # init of reaction indexes form trigable reactions
        r_ind=[]
        
        # generates of reaction indexes form trigable reactions
        for i in range(self.MpDf.shape[1]):
             if (sp & self.ReacListBt[i]) == self.ReacListBt[i]:
                 r_ind.append(i)
        
        if len(r_ind)==0:
            return True
        
        # init of  indexes form trigable reactions
        sp_ind=[]
        for i in range(len(sp)):
            if sp[i]==1:
                sp_ind.append(i)
        
        # relative Soichiometric matrix to present species and reaction
        S=np.array(self.MrDf-self.MpDf)
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
        
        if ret_vect:
            
            v=np.zeros(self.MpDf.shape[1])
            v[r_ind]=res.x
            return res.success , v
        
        return res.success
        
    
    # Subset funcntion for bitarrays
    def BtIsSubsetofBt(self, a, b):
        return (a & b) == a
      
    # Function that returns Bitarray indexes for postition with 1 value 
    def getIndArrayFromBt(self,bt):
        ind=[]
        for i in range(len(bt)):
            if bt[i]==1:
                ind.append(i)
        
        return ind
    
    # Function that returns Bitarray for postition from a 0 and 1 array 
    def getBtFromIndArray(self,ind,size):
        btarray=bt(size)
        btarray.setall(0)
        for i in ind:
           btarray[i]=1
        
        return btarray
    
    # Function that receives a set of species (sp_set) and returns 
    # the idexes of which are outflow.
    def getOutflowFromSp(self,sp_set,set_type=False):
                
        # Checks if input is or not bitarray, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
        
        inf=sp.copy()
        inf.setall(0)
        for i in range(len(self.ReacListBt)):
            if (self.ReacListBt[i].any()) and (not self.ProdListBt[i].any()):
                inf |= self.ReacListBt[i]
        
        inf&=sp
        
        if set_type:
            return self.SpIdStrArray[self.getIndArrayFromBt(inf)]
        
        return self.getIndArrayFromBt(inf)
        
    # Function that receives a set of species (sp_set) and returns 
    # the idexes of which are inflow.
    def getInflowFromSp(self,sp_set,return_type="index"):
                
        # Checks if input is or not bitarray, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
        
        out=sp.copy()
        out.setall(0)
        for i in range(self.MpDf.shape[1]):
            if (self.ProdListBt[i].any()) and (not self.ReacListBt[i].any()):
                out |= self.ProdListBt[i]
        
        out&=sp
        if return_type=="id":
            return self.SpIdStrArray[self.getIndArrayFromBt(out)]
        
        elif return_type=="index":
            return self.getIndArrayFromBt(out)
        
        elif return_type=="bt":
            return out
    
    def getRandNat(self,Sum,size):
        '''
        

        Parameters
        ----------
        Sum : int > 0
            How much the components sums.
        size : int > 0
            size of the vector.

        Returns
        -------
        A vector of naturals of size size, hows componennts sums Sum.

        '''
        
        if isinstance(Sum, int) and isinstance(size, int) and (size>0) and (Sum>0):
            return np.random.multinomial(Sum, [1/float(size)] * size)
                    
        else:
            raise TypeError("Only naturals numbers are allowed")
            
        
    # random generation of reaction networks
    @classmethod
    def setSimpleRandomgenerate(cls,Ns=12,rv=None):
        '''
        
    
        Parameters
        ----------
        Ns : int, optional
            Number of species. The default is 12.
        rv : vector or int, optional
            Either a natural number or a vector of dimension 7. 
            If rv is a natural number, it will
            be used to generate a randomized vector of dimension 7 for the total number
            of reactions rv
            
            If rv is a vector of
            dimension 7, the network will contain rv[i] reactions of type i, where type i is:
            i=0. inflow ∅ → x
            i=1. outflow x → ∅
            i=2. transformation x → y
            i=3. synthesis x + y → z
            i=4. decomposition z → x + y
            i=5. single replacement x + y → x + z
            i=6. double replacement x + y → z + w
           
           Finally if rv =None, rv is a random vector of dimension 7 and
            total components sum equals to Ns. The default is None.
    
        Returns
        -------
        Simple generated random network.
    
        '''
        
        # Chacking consistency of number of species
        if (Ns > 0) and isinstance(Ns, int):
            
            if rv is None :
                rv=np.random.multinomial(Ns, [1/7.] * 7)
                
            elif isinstance(rv,int) and (rv>0):
                rv=np.random.multinomial(rv, [1/7.] * 7)
            elif (len(rv)<7) or (not all(isinstance(n, int) for n in rv)):
                raise TypeError("inconsistency of input variables")
        else:
            raise TypeError("Ns in not a natural number")
        

        # Creating species 
        sp=list(map(lambda x: "s"+str(x+1).zfill(1+int(np.floor(np.log10(Ns)))),range(Ns)))

        #Cratinng
        Nr=sum(rv)

        
        mr=np.zeros([Nr,Ns])
        mr=pd.DataFrame(mr,columns=sp)

        mp=copy.copy(mr)
        
        k=0 #reaction counter
        # Adding type 0 reactions
        if rv[0]!=0:

            rsp=np.random.choice(Ns, size=rv[0], replace=False)
            
            for i in rsp:
                mp.iloc[k,i]=1
                k+=1
    
        # Adding type 1 reactions
        if rv[1]!=0:

            rsp=np.random.choice(Ns, size=rv[1], replace=False)
            
            for i in rsp:
                mr.iloc[k,i]=1
                k+=1
            
        # Adding type 2 reactions
        if rv[2]!=0:

            l=0
            while l<rv[2]:
            
                rsp=np.random.choice(Ns, size=2, replace=False)
                
                reac=np.zeros(Ns)
                reac[rsp[0]]=1
                prod=np.zeros(Ns)
                prod[rsp[1]]=1
                if not ((mr == reac).all(1).any() and (mp == prod).all(1).any()):
                    
                    mr.iloc[k,rsp[0]]=1
                    mp.iloc[k,rsp[1]]=1
                    l+=1
                    k+=1
        
        # Adding type 3 reactions
        if rv[3]!=0:

            l=0
            while l<rv[3]:
            
                rsp=np.random.choice(Ns, size=3, replace=False)
                
                reac=np.zeros(Ns)
                reac[rsp[0]]=1
                reac[rsp[1]]=1
                prod=np.zeros(Ns)
                prod[rsp[2]]=1
                if not ((mr == reac).all(1).any() and (mp == prod).all(1).any()):
                    
                    mr.iloc[k,rsp[0]]=1
                    mr.iloc[k,rsp[1]]=1
                    mp.iloc[k,rsp[2]]=1
                    l+=1
                    k+=1
        
        # Adding type 4 reactions
        if rv[4]!=0:
            
            l=0
            while l<rv[4]:
            
                rsp=np.random.choice(Ns, size=3, replace=False)
                
                reac=np.zeros(Ns)
                reac[rsp[0]]=1
                prod=np.zeros(Ns)
                prod[rsp[1]]=1
                prod[rsp[2]]=1
                if not ((mr == reac).all(1).any() and (mp == prod).all(1).any()):
                    
                    mr.iloc[k,rsp[0]]=1
                    mp.iloc[k,rsp[1]]=1
                    mp.iloc[k,rsp[2]]=1
                    l+=1
                    k+=1

        # Adding type 5 reactions
        if rv[5]!=0:
            
            l=0
            while l<rv[5]:
            
                rsp=np.random.choice(Ns, size=3, replace=False)
                
                reac=np.zeros(Ns)
                reac[rsp[0]]=1
                reac[rsp[1]]=1
                prod=np.zeros(Ns)
                prod[rsp[0]]=1
                prod[rsp[2]]=1
                if not ((mr == reac).all(1).any() and (mp == prod).all(1).any()):
                    
                    mr.iloc[k,rsp[0]]=1
                    mr.iloc[k,rsp[1]]=1
                    mp.iloc[k,rsp[0]]=1
                    mp.iloc[k,rsp[2]]=1
                    l+=1
                    k+=1
                    
        # Adding type 6 reactions
        if rv[6]!=0:
            
            l=0
            while l<rv[6]:
            
                rsp=np.random.choice(Ns, size=4, replace=False)
                
                reac=np.zeros(Ns)
                reac[rsp[0]]=1
                reac[rsp[1]]=1
                prod=np.zeros(Ns)
                prod[rsp[2]]=1
                prod[rsp[3]]=1
                if not ((mr == reac).all(1).any() and (mp == prod).all(1).any()):
                    
                    mr.iloc[k,rsp[0]]=1
                    mr.iloc[k,rsp[1]]=1
                    mp.iloc[k,rsp[2]]=1
                    mp.iloc[k,rsp[3]]=1
                    l+=1
                    k+=1
                    
            
        mr=mr.T
        mp=mp.T
        
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
            
        #sorting mp and mr dataframes
        sp_i=mp.index.values.tolist()
        sorted_sp_i=mp.sort_index(axis=0).index.values.tolist()
        sorted_ind = [sp_i.index(i) for i in sorted_sp_i]
        
        for i in range(len(prod)):
            tmp_prod=bt(len(mr.index))
            tmp_reac=bt(len(mr.index))
            for j in range(len(mr.index)):
                tmp_prod[j]=prod[i][sorted_ind[j]]
                tmp_reac[j]=reac[i][sorted_ind[j]]
            
            prod[i]=tmp_prod
            reac[i]=tmp_reac
        
        mp=mp.sort_index(axis=0)
        mr=mr.sort_index(axis=0)
        
        out =cls()
        out.MrDf=mr
        out.MpDf=mp
        out.SpIdStrArray=np.array(mr.index)
        out.SpNameStrArray=out.SpIdStrArray.copy()
        out.ReacListBt=reac
        out.ProdListBt=prod
        out.FilenameStr=None
        out.IsTextBool=False
        out.IsSbmlBool=False
        
        return out
            
    # the most simple random network generator, Nr reactions (>1), Ns species (>1)
    # a minimal reaction network is randomly created where each reaction has one reactant and one product and each species
    # is used at least once as reactant and once as product, extra assignments are carried out randomly over reactions
    # there are no inflow or outflow reactions, redundant or null reactions may be generated (rn.merge will filter them out)
    # dist is a log scaled distribution in the [-1,1] range representing locality
    # pr and pp are a log scaled penalization for the repeated use of species as reactants or products
    @classmethod
    def setRandomgeneratedNoInflow(cls,Nr=12,Ns=None,extra=.4, dist=lambda x: x*0+1, pr=0, pp=None):
        
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
        sp=list(map(lambda x: "s"+str(x+1).zfill(1+int(np.floor(np.log10(Ns)))),range(Ns)))
        mr=np.zeros([Nr,Ns])
        mr=pd.DataFrame(mr,columns=sp)
        mr=mr.T
        mp=copy.copy(mr)
        # (2) assignment of one random reactant and one random product (different from the reactant) to each reaction
        for i in np.random.choice(range(Nr),Nr,False):
            d=np.array(list(map(dist,list(map(lambda x: rnorm(x-xr[i]),xs)))))
            dr = d - pr*usr; 
            dr = np.exp(dr-np.max(dr))
            sr = np.random.choice(range(Ns),1,p=norm(dr))
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
          r = np.random.choice(range(Nr),1,p=norm(np.exp(d)))
          mr.iloc[s,r] = mr.iloc[s,r] + 1
        
        # (4) assignment of species not used as products to a random reaction (eventually the same reaction)
        i = np.where(usp==0)[0]
        
        for s in np.random.choice(i,len(i)):
          d=np.array(list(map(dist,list(map(lambda x: rnorm(x-xs[s]),xr)))))
          r = np.random.choice(range(Nr),1,p=norm(np.exp(d)))
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
            
        #sorting mp and mr dataframes
        sp_i=mp.index.values.tolist()
        sorted_sp_i=mp.sort_index(axis=0).index.values.tolist()
        sorted_ind = [sp_i.index(i) for i in sorted_sp_i]
        
        for i in range(len(prod)):
            tmp_prod=bt(len(mr.index))
            tmp_reac=bt(len(mr.index))
            for j in range(len(mr.index)):
                tmp_prod[j]=prod[i][sorted_ind[j]]
                tmp_reac[j]=reac[i][sorted_ind[j]]
            
            prod[i]=tmp_prod
            reac[i]=tmp_reac
        
        mp=mp.sort_index(axis=0)
        mr=mr.sort_index(axis=0)
        
        out =cls()
        out.MrDf=mr
        out.MpDf=mp
        out.SpIdStrArray=np.array(mr.index)
        out.SpNameStrArray=out.SpIdStrArray.copy()
        out.ReacListBt=reac
        out.ProdListBt=prod
        out.FilenameStr=None
        out.IsTextBool=False
        out.IsSbmlBool=False
        
        return out

    
    def setExtraRandomgenerated(self,p=.1,m=2,Nse=None,extra=None,l="x"):
        
        mr = self.MrDf
        mp = self.MpDf
        Ns = mr.shape[0]
        Nr = mr.shape[1]
        
        if Nse is None:
            Nse=np.ceil(Ns*p).astype(int)
        if extra is None:
            extra=np.round(m*Nr).astype(int)

        # (1) adding extra species
        me=np.zeros((Nr,Nse))
        
        # fidin if prefix l is alredy used.
        add_species=self.SpIdStrArray[list(map(lambda x: x.find(l)!=-1,self.SpIdStrArray))]
        if any(add_species):
            begin_idx=max(list(map(lambda x: int(x.replace("x","")),add_species)))
        else:
            begin_idx=0
        
        me=pd.DataFrame(me,columns=list(map(lambda x: l+str(x+1).zfill(1+int(np.floor(np.log10(Nse)))),range(begin_idx,Nse+begin_idx))))
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
        self.MrDf=mr
        self.MpDf=mp
        self.SpIdStrArray=np.array(mr.index)
        self.SpNameStrArray=self.SpIdStrArray.copy()
        self.ReacListBt=reac
        self.ProdListBt=prod
    
    # function that adds a percentage of additional (extra) inflow 
    # reactions to an existing network 
    def setExtraRandomgeneratedInflow(self,extra=0.1):
        Ns = self.MrDf.shape[0]
        
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
        mr.index=self.MrDf.index.copy()
        mp=pd.DataFrame(mp)
        mp.index=self.MrDf.index.copy()
        
        # adding reaction to respective patices
        self.MrDf = pd.concat([self.MrDf,mr],axis=1)
        self.MrDf.columns=range(self.MrDf.shape[1])
        self.MpDf = pd.concat([self.MpDf,mp],axis=1)
        self.MpDf.columns=range(self.MpDf.shape[1])
        
        # creating extra bitset variables for the support and products
        reac=[]
        prod=[]
        for i in range(self.MpDf.shape[1]):
            
            r_sp=bt(self.MrDf.shape[0])
            r_sp.setall(0)
            p_sp=r_sp.copy()
            
            for j in np.where(self.MrDf.iloc[:,i]!=0)[0]:
                r_sp[j]=1
            for j in np.where(self.MpDf.iloc[:,i]!=0)[0]:
                p_sp[j]=1
            
            reac.append(r_sp)
            prod.append(p_sp)       
        
        self.ReacListBt=reac
        self.ProdListBt=prod
    
    def addInflow(self,sp_set):
        '''
        

        Parameters
        ----------
        sp_set : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        
        # selection of species that will considere as inflow species
        SpIds = sp_set.search(1)
        inflow_sp = self.getInflowFromSp(self.SpIdStrArray)
        
        SpIds =list(set(SpIds)-set(inflow_sp))
        
        
        
        k=self.MpDf.shape[1]
        for l in SpIds:
            self.MpDf[k]=0.0
            self.MrDf[k]=0.0
            self.MpDf.loc[self.SpIdStrArray[l],k]=1.0
            
            p_bt=bt(self.MpDf.shape[0])
            p_bt.setall(0)
            r_bt=p_bt.copy()
            
            self.ReacListBt.append(r_bt)
            
            p_bt[l]=1
            self.ProdListBt.append(p_bt)
            k+=1
    
    # function that adds a percentage of additional (extra) inflow 
    # reactions to an existing network
    def setExtraRandomgeneratedOutflow(self,extra=0.1):
        Ns = self.MrDf.shape[0]
        
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
        mr.index=self.MrDf.index.copy()
        mp=pd.DataFrame(mp)
        mp.index=self.MrDf.index.copy()
        
        # adding reaction to respective patices
        self.MrDf = pd.concat([self.MrDf,mr],axis=1)
        self.MrDf.columns=range(self.MrDf.shape[1])
        self.MpDf = pd.concat([self.MpDf,mp],axis=1)
        self.MpDf.columns=range(self.MpDf.shape[1])
        
        # creating extra bitset variables for the support and products
        reac=[]
        prod=[]
        for i in range(self.MpDf.shape[1]):
            
            r_sp=bt(self.MrDf.shape[0])
            r_sp.setall(0)
            p_sp=r_sp.copy()
            
            for j in np.where(self.MrDf.iloc[:,i]!=0)[0]:
                r_sp[j]=1
            for j in np.where(self.MpDf.iloc[:,i]!=0)[0]:
                p_sp[j]=1
            
            reac.append(r_sp)
            prod.append(p_sp)       
        
        self.ReacListBt=reac
        self.ProdListBt=prod

    
    # A wrapper of function setRandomgeneratedNoInflow, but adding percentage of input inflow 
    # and input outflow reactions.
    @classmethod
    def setRandomgeneratedWithInflow(cls,Nr=12,Ns=None,extra=.4, dist=lambda x: x*0+1, pr=0, pp=None, inflow=0.1, outflow=0.1):
        
        # inizialization of 
        out=cls.setRandomgeneratedNoInflow(Nr,Ns,extra, dist, pr, pp)
        out.setExtraRandomgeneratedInflow(inflow)
        out.setExtraRandomgeneratedOutflow(outflow)
        out.SetRnClean()
        
        return out
            
    # Function that cleans up reaccion redundancies and unsed species.
    def SetRnClean(self):
        
        mr=self.MrDf.copy()
        mp=self.MpDf.copy()
        
        i = np.where(mr.sum(axis=1)+mr.sum(axis=1)==0)[0]
        if len(i)>0:
            mr = mr.drop(mr.index[i],axis=0,inplace=False)  # unused species are eliminated
            mp = mp.drop(mp.index[i],axis=0,inplace=False)  # unused species are eliminated
                
        i = np.where(list(map(lambda k: all(mr.iloc[:,k]==mp.iloc[:,k]), range(mr.shape[1]))))[0]

        if len(i)>0:
            mr = mr.drop(mr.columns[i],axis=1,inplace=False)  # unused species are eliminated
            mp = mp.drop(mp.columns[i],axis=1,inplace=False)  # unused species are eliminated
          
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
        self.MrDf=mr
        self.MpDf=mp
        self.SpIdStrArray=np.array(mr.index)
        self.SpNameStrArray=self.SpIdStrArray.copy()
        self.ReacListBt=reac
        self.ProdListBt=prod
    
    @classmethod    
    def setRandomgeneratedBoolean(cls,init_Nr=10,add_Nr=10,bsp_N=10,bsp=None,bsp_wf=1/3,bsp_w=None,
                                  sp_len_max=10, sp_paired_p=None,
                                  sp2_wf=1/3, sp2_len_w=None, sp2_right_side_p=.5,
                                  new_sp_len_w=None,synt_reac_p=.5,close=True,reuse=False):
    
        
        RN=gen_rn_csp(bsp_N, bsp ,bsp_wf ,bsp_w, sp_len_max, sp_paired_p, 
                      sp2_wf, sp2_len_w, sp2_right_side_p, new_sp_len_w, synt_reac_p)
        out=cls()
        RN.test(init_Nr,add_Nr,close,reuse)
        #RN.add_reactions(Nr,close,reuse)
        mr, mp =RN.get_matrices()
        out.MpDf=pd.DataFrame(mr,index=np.array(RN.sp))
        out.MrDf=pd.DataFrame(mp,index=np.array(RN.sp))
        out.SpIdStrArray=np.array(RN.sp)
        out.SpNameStrArray=out.SpIdStrArray.copy()

        reac=[]
        prod=[]
        for i in range(out.MpDf.shape[1]):
            
            r_sp=bt(out.MpDf.shape[0])
            r_sp.setall(0)
            p_sp=r_sp.copy()
            
            for j in np.where(out.MrDf.iloc[:,i]!=0)[0]:
                r_sp[j]=1
            for j in np.where(out.MpDf.iloc[:,i]!=0)[0]:
                p_sp[j]=1
            
            reac.append(r_sp)
            prod.append(p_sp)
            
        out.ReacListBt=reac
        out.ProdListBt=prod
        out.FilenameStr=None
        out.IsTextBool=False
        out.IsSbmlBool=False
        
        return out
        
    
 