#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 20:12:37 2022

@author: pmaldona

Closed Reaction Network Strucutre Library
"""
from .RNIRG import RNIRG
import numpy as np
import pandas as pd
from bitarray import bitarray as bt
from bitarray import frozenbitarray as fbt
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
from bitarray.util import subset
from collections import Counter as count
import json
import pickle as pk
from pulp import *
import re
from functools import reduce

# Class that calculates the synergistic and organizational structure 
# of the closed reactants of a reaction network. 
class CRNS(RNIRG):
    
    # Function that generates the basic sets of a reaction network. When 
    # executed, it generates the following members of the class, which 
    # correspond to a list of bitarrays, where each element of the list 
    # corresponds to a basic set or generator as the case may be:
    
    # BSpListBt: species contained in each basic molecule set
    # BRpListBt: reactions contained in each basic molecule set
    # GSpListBt: species contained in each generator
    # GRpListBt: reactions contained in each generator 
    # BReacSpListBt: species of reactants contained in each basic molecule set 
    # BProdSpListBt: species of products contained in each basic molecule set
    # BStoichioPositiveSpListBt: positive stoichiometric species contained in basic molecule set
    # BStoichioNegativeSpListBt: negative stoichiometric species contained in basic molecule set
    # GInBListBt: generator contained in each basic molecule set
    # ConnectedBListBt: basic molecule sets that are connected (whose has a non empty intersection 
    # between the support and the products o viceversa).
    # NotContainedBListBt: basic molecule sets that are not contained (whose has a non empty intersection 
    # between the support and the products o viceversa).      
    # RpInWhichGArray: corresponds to a array of integers, of the length 
    # of the number of reactions, where the integer indicates to which 
    # generator the reaction belongs.
    def setGenerators(self,verbose=False):
        # creating all closed set for closure of reactant part
        c_reac=[]
          
        k=1
        for i in self.ReacListBt:
            if verbose:
                print(k," closures of ",len(self.ReacListBt))
            c_reac.append(self.getClosureFromSp(i,bt_type=True))
            k+=1
            
        # number of equivalence clases
        xeqc=0
        # class equivalence vector
        x_r_a=-np.ones(self.MpDf.shape[1])
        x_r_a=x_r_a.tolist()
        # list of basic sets
        sp_b=[]
        
        # number of steps
        st=0
        # creating equivalence classes for each reaction that generate the same
        # closure and so each basic can be created
        
        for i in range(len(x_r_a)):
            if verbose:
                print(i+1," equivalance classes of ",len(x_r_a))
            
            st+=1
            # reaction already assigned to equivalnce class
            if x_r_a[i]>=0:
               continue
           
            x_r_a[i]=xeqc
            sp_b.append(c_reac[i])
            
            for j in range(i,self.MpDf.shape[1]):
                st+=1
                if c_reac[i]==c_reac[j]:
                    x_r_a[j]=xeqc
                    
            xeqc+=1
        
        self.BSpListBt=sp_b
        self.RpInWhichGArray=x_r_a
        
        # assingnig species contained each generator
        sp_a=[]
        
        b_sp=bt(self.MpDf.shape[0])
        b_sp.setall(0)
        
        for i in range(xeqc):
            sp_a.append(b_sp.copy())
            
        for i in range(self.MpDf.shape[1]):
            sp_a[x_r_a[i]]|= self.ReacListBt[i]|self.ProdListBt[i] 
        
        
        self.GSpListBt=sp_a
        
        
        # assigning reaction contained each generator
        r_a=[]
        b_r=bt(self.MpDf.shape[1])
        b_r.setall(0)
        
        for i in range(xeqc):
            for j in np.where(np.array(x_r_a)==i)[0]:
                b_r[j]=1   
            
            r_a.append(b_r.copy())
            b_r.setall(0)
            
            
        self.GRpListBt=r_a
        
        # Creation of other related variables:
        a_b=[] # generators (equivalence classes) contained in each closure
        r_b=[] # reactions supported by each basic (equivalence class)
        rsp_b=[] # reactants species contained each basic (equivalence class)
        psp_b=[] # product species contained each basic (equivalence class)
        sp_sp_b=[] # stoichiometric positive species contained in the closure of each generator
        sn_sp_b=[] # stoichiometric negative species contained in the closure of each generator
        
        # Assignation of other related variables:
        for i in range(len(self.BSpListBt)):
            if verbose:
                print(i+1," bitarrays assignations of ",len(self.BSpListBt))
            a_b.append(self.getGBtInSpBt(self.BSpListBt[i]))
            r_b.append(b_r.copy())
            rsp_b.append(b_sp.copy())
            psp_b.append(b_sp.copy())
            sp_sp_b.append(b_sp.copy())
            sn_sp_b.append(b_sp.copy())
            
            for j in self.getIndArrayFromBt(a_b[i]):
                r_b[i]|=r_a[j]
                for k in self.getIndArrayFromBt(r_a[j]): 
                    rsp_b[i]|=self.ReacListBt[k]
                    psp_b[i]|=self.ProdListBt[k]
                    sp_sp_b[i]|=bt(self.MrDf.iloc[:,k] < self.MpDf.iloc[:,k])
                    sn_sp_b[i]|=bt(self.MrDf.iloc[:,k] > self.MpDf.iloc[:,k])
            
        self.GInBListBt=a_b
        self.BRpListBt=r_b
        self.BReacSpListBt=rsp_b
        self.BProdSpListBt=psp_b
        self.BStoichioPositiveSpListBt=sp_sp_b
        self.BStoichioNegativeSpListBt=sn_sp_b
       
        # Creating connecting basics and dynamically connected basics sets
        conn=[]
        dyn_conn=[]
        
        # Creating connecting basics and dynamically connected basics sets
        conn=[]
        dyn_conn=[]
        
        b_c=bt(len(self.BSpListBt))
        b_c.setall(0)
        
        for i in range(len(self.BSpListBt)):
            
            conn.append(b_c.copy())
            dyn_conn.append(b_c.copy())
            
        for i in range(len(conn)):
            if verbose:
                print(i," connectivity assingments of ",len(self.BSpListBt))
            for j in range(len(conn)):
                # connectivity conditions
                if (not j==i) and a_b[i][j]==0 and (rsp_b[i] & psp_b[j]).any():
                    dyn_conn[i][j]=1
                    dyn_conn[j][i]=1
                elif (not j==i) and a_b[i][j]==0 and (psp_b[i] & rsp_b[j]).any():
                    dyn_conn[i][j]=1
                    dyn_conn[j][i]=1
                # Hasse condition (all connected)
                if (not j==i) and a_b[i][j]==0:
                    conn[i][j]=1
                    # conn[j][i]=1
        
        self.NotContainedBListBt=conn
        self.ConnectedBListBt=dyn_conn
        return(st)
    
            
    # Function that returns the number of basic sets in which each species of 
    # the vector sp appears. The input vector sp can be a list of 
    # strings or a bitarray.
    def getSpPresenceInBGArray(self,sp_set):
        # Checks if input is or not bitarray, if it's no, it make the 
        # transformation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
        
        sp_b_presc=np.zeros(len(sp))
        sp_a_presc=np.zeros(len(sp))
        for i in sp.itersearch(1):
            for j in range(len(self.BSpListBt)):
                if self.BSpListBt[j][i]==1:
                    sp_b_presc[i]+=1
            for j in range(len(self.GSpListBt)):
                if self.GSpListBt[j][i]==1:
                    sp_a_presc[i]+=1
        sp_b_presc=sp_b_presc[self.getIndArrayFromBt(sp)]
        sp_a_presc=sp_a_presc[self.getIndArrayFromBt(sp)]
        
        sp_b_presc=pd.DataFrame(sp_b_presc)
        sp_b_presc=sp_b_presc.transpose()
        sp_b_presc.columns=self.SpIdStrArray[self.getIndArrayFromBt(sp)]
        
        sp_a_presc=pd.DataFrame(sp_a_presc)
        sp_a_presc=sp_a_presc.transpose()
        sp_a_presc.columns=self.SpIdStrArray[self.getIndArrayFromBt(sp)]
        
        return([sp_b_presc,sp_a_presc])
    
    
    # Function that returns the number of basic sets in which each reaction of 
    # the vector r appears. The input vector v can be a list of 
    # strings or a bitarray.
    def getRpPresenceInBArray(self,r_set):
        # Checks if input is or not bitarray, if it's no, it make the 
        # transformation
        if not (isinstance(r_set,bt)):
            v=bt(self.MpDf.shape[1])
            v.setall(0)
            
            for i in r_set:
                v[i]=1
        else:
            v=r_set

        
        r_b_presc=np.zeros(len(v))
        
        for i in v.itersearch(1):
            for j in range(len(self.BRpListBt)):
                if self.BRpListBt[j][i]==1:
                    r_b_presc[i]+=1

        r_b_presc=r_b_presc[self.getIndArrayFromBt(v)]
        
        r_b_presc=pd.DataFrame(r_b_presc)
        r_b_presc=r_b_presc.transpose()
        
        names=[]
        
        for i in self.MpDf.columns:
            names.append("r"+str(i))
        
        r_b_presc.columns=np.array(names)[self.getIndArrayFromBt(v)]
          
        return(r_b_presc)
    
    # Function that returns the number of basic sets in which each species of 
    # the vector sp appears. The input vector sp can be a list of 
    # strings or a bitarray.
    def plotSpPresenceInBG(self,sp_set):
        histo=self.getSpPresenceInBGArray(sp_set)
        
        names=histo[0].columns.to_list()
        values=histo[0].iloc[0].tolist()
    
        fig = plt.figure(figsize = (10, 5))
 
        # creating the bar plot
        plt.bar(names, values, color ='maroon',
                width = 0.4)
         
        plt.xlabel("Species")
        plt.ylabel("Number of basics sets")
        plt.title("Number of basic sets that contain each species")
        plt.show()
        
        histo=self.getSpPresenceInBGArray(sp_set)
        
        names=histo[1].columns.to_list()
        values=histo[1].iloc[0].tolist()
    
        fig = plt.figure(figsize = (10, 5))
 
        # creating the bar plot
        plt.bar(names, values, color ='maroon',
                width = 0.4)
         
        plt.xlabel("Species")
        plt.ylabel("Number of generators")
        plt.title("Number of generators sets that contain each species")
        plt.show()
    
    
    # Function that plots the number of basic sets in which each reaction of 
    # the vector r appears. The input vector v can be a list of 
    # strings or a bitarray.
    def plotRpPresenceInB(self,v):
        
        histo=self.getRpPresenceInBArray(v)
        
        names=histo.columns.to_list()
        values=histo.iloc[0].tolist()
    
        fig = plt.figure(figsize = (10, 5))
 
        # creating the bar plot
        plt.bar(names, values, color ='maroon',
                width = 0.4)
         
        plt.xlabel("Reactions")
        plt.ylabel("Number of Basics")
        plt.title("Number of Basic that contain each reaction")
        plt.show()
    
    
    # Function of species that returns the generator that contains the species
    def getGBtInSpBt(self, sp):
        p=bt(len(self.BSpListBt))
        p.setall(0)
        for i in range(len(self.GSpListBt)):
            if (sp & self.GSpListBt[i]) == self.GSpListBt[i]:
                p[i]=1
        
        return p
    
    
    # Function that returns the contains the species contained in a Bitarray ofgetGBtContribBBt
    # generators
    def getSpBtInGBt(self, p):
        sp=bt(len(self.SpIdStrArray))
        sp.setall(0)
        for i in self.getIndArrayFromBt(p):
            
                sp|=self.BSpListBt[i]
        
        return sp
    
    
    # For a given bitarray of contained basic set, return the connected basics 
    # (basics that are not contained)
    def getGBtNotInBBt(self,s):
        
        c=bt(len(self.BSpListBt))
        c.setall(0)
        
        for i in self.getIndArrayFromBt(s):
            c|=self.NotContainedBListBt[i]
            
        c= c & ~s
        return c
        
    
    # For a given bitarray of contained basic set, return the dynamically 
    # connected basics sets to it
    def getGBtConnectedToBBt(self,s,inculde_itself=False):
        c=bt(len(self.BSpListBt))
        c.setall(0)
        
        for i in self.getIndArrayFromBt(s):
            c|=self.ConnectedBListBt[i]
            
        if not inculde_itself:
            c= c & ~s
        else:
            c= c | s
        return c
    
    def setTransGConn(self):
        
        trans_dyn_conn=[]
        dyn_conn_bt=self.GInBListBt[0].copy()        
        for i in self.ConnectedBListBt:
            dyn_conn_bt=self.getGBtConnectedToBBt(i,inculde_itself=True)
            closure=False
            while closure:
                add_conn_bt=self.getGBtConnectedToBBt(dyn_conn_bt,inculde_itself=True)
                if subset(add_conn_bt,dyn_conn_bt):
                    closure=True
                dyn_conn_bt=add_conn_bt
                
            trans_dyn_conn.append(dyn_conn_bt)
            
        self.TransConnectedBListBt=trans_dyn_conn
            
    # For a given bitarray of contained basic set, return the dynamically 
    # connected basics sets to it
    def getGBtTransConnectedToBBt(self,s,inculde_itself=False):
        c=bt(len(self.BSpListBt))
        c.setall(0)
        
        for i in self.getIndArrayFromBt(s):
            c|=self.TransConnectedBListBt[i]
            
        if not inculde_itself:
            c= c & ~s
        else:
            c= c | s
        return c
    
    # For a given bitarray of contained basic set, return the basic sets that 
    # can contribute for the current set to be semi-self-maintained 
    def getGBtContribBBt(self,s):
        
        # negative sotichiometric species
        
        n=bt(self.MpDf.shape[0])
        n.setall(0)
        p=n.copy()
        pp=n.copy()
        
        for i in self.getIndArrayFromBt(s):
            p|=self.BStoichioPositiveSpListBt[i]
            n|=self.BStoichioNegativeSpListBt[i]
                       
        # p correspond to bitarray whit positive stoichiometry
        n= n & ~p
        
        c=s.copy()
        if not n.any():
           c.setall(0)
           return c
       
        # possible sets that can contribute to be semi-self-maintained
        c = ~c
        
        # eliminating basics form c that no will contribute
        for i in self.getIndArrayFromBt(c):
            if not (self.BStoichioPositiveSpListBt[i] & n).any():
                c[i]=0
            else:
                pp|=self.BStoichioPositiveSpListBt[i]
        
        # pp are possible produces species if they do not fulfill consumed 
        # species, then there is no contribution to be a semi-self-maintained

        if (n & ~pp).any():
            c.setall(0)
            return c
        else:
            return c
        
    
    # Synergistic structure calculation function, requires the setGenerators() 
    # function to be executed beforehand .It returns an directed multi-graph
    # type networkx (SynStrNx), where each node correspond to hashed bitarray of 
    # contained basic basic sets. Each node also contain attributes such as level (level), 
    # set of contained species (sp) and if the set is semi-self-maintained (ssm)
    # and an organization (is_org). 
    # The edges characterize the closures with the different basic sets, where 
    # the arrival set corresponds to the closure, the key to the 
    # basic with which the closure was performed and the attribute whether the closure 
    # is synergistic or not. 
    # The function also creates the class member (SynStrSsmListSpArray) and (SynStrOrgListSpArray) which 
    # corresponds to a list of the closed semi-self-maintained set and 
    # organizations respectively
    def setSynStr(self,partial_save=None,verbose=False):
        if not hasattr(self, 'GInBListBt'):
            print("The basic sets have not been initialized, please run the setGenerators() function.")
            return 
        # Initialization of the synergistic structure as a multigraph object
        G = nx.MultiDiGraph()
        # Initialization of the list of semi-self-maintained sets and organizations
        ssms=[]
        ssms_bt=[]
        org=[]
        org_bt=[]
        # step measure
        st=0
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.GInBListBt)):
            st+=1
            is_ssm=self.isSsmFromSp(self.BSpListBt[i])
            if is_ssm:
                is_org=self.isSmFromSp(self.BSpListBt[i])
                G.add_node(fbt(self.GInBListBt[i]),level=self.GInBListBt[i].count(),
                       sp=self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])],
                       is_ssm=is_ssm,is_org=is_org,is_basic=True,basic_id=i)
                ssms.append(self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])])
                ssms_bt.append(self.BSpListBt[i])
                if is_org:
                    org.append(self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])])
                    org_bt.append(self.BSpListBt[i])
            else:           
                G.add_node(fbt(self.GInBListBt[i]),level=self.GInBListBt[i].count(),
                       sp=self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])],
                       is_ssm=False,is_org=False,is_basic=True,basic_id=i)

                
        # Generation of the multigraph of the synergistic structure by set level (contained basics)
        for i in range(len(self.BSpListBt)):
            
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
            n_nodes=len(nodes)
            node_n=0
             # Generating closures whit connected basics sets for each set in level i
            for j in nodes:
                node_n+=1
                if verbose:
                    print("level: ",i+1, "from ",len(self.BSpListBt),", node: ",node_n," from ",n_nodes)
                for k in self.getIndArrayFromBt(self.getGBtNotInBBt(bt(j))):
                    st+=1    
                    # Closure result
                    cr_sp=self.getClosureFromSp(self.getSpBtInGBt(bt(j) | self.GInBListBt[k]),bt_type=True)
                    cr_a=fbt(self.getGBtInSpBt(cr_sp))
                    
                    # node is added if is not in structrue
                    if not (cr_a in G):
                        is_ssm=self.isSsmFromSp(cr_sp)
                        if is_ssm:
                            is_org=self.isSmFromSp(cr_sp)
                            G.add_node(cr_a,level=cr_a.count(),
                                   sp=self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)],
                                   is_ssm=is_ssm,is_org=is_org,is_basic=False)
                            ssms.append(self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)])
                            ssms_bt.append(cr_sp)
                            if is_org:
                                org.append(self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)])
                                org_bt.append(cr_sp)
                        else:           
                            G.add_node(cr_a,level=cr_a.count(),
                                   sp=self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)],
                                   is_ssm=False,is_org=False,is_basic=False)
                        
                    # Adding edges corresponding to the colsure, and verifing if is a sinergy:
                    if cr_a.count() > (bt(j)|self.GInBListBt[k]).count():
                       G.add_edge(j,cr_a,key=fbt(self.GInBListBt[k]),syn=True,added_basic=k)
                    else:
                       G.add_edge(j,cr_a,key=fbt(self.GInBListBt[k]),syn=False,added_basic=k)
                   # if cr_a.count() > (bt(j)|self.GInBListBt[k]).count():
                   #    G.add_edge(j,cr_a,key=fbt(self.GInBListBt[k]),syn=True,added_basic=k)
                   # else:
                   #    G.add_edge(j,cr_a,key=fbt(self.GInBListBt[k]),syn=False,added_basic=k)
           
            if not partial_save is None:
                self.saveToPkl(partial_save)
        
        self.SynStrNx=G
        self.SynStrSsmListSpArray=ssms
        self.SynStrOrgListSpArray=org
        self.SynStrSsmListBtArray=ssms_bt
        self.SynStrOrgListBtArray=org_bt
        return(st)
     
    def getconnStrNx(G):
        H=nx.MultiDiGraph()
        
        
                

    # Synergistic structure calculation function, requires the setGenerators() 
    # function to be executed beforehand .It returns an directed multi-graph
    # type networkx (setSsmStr), where each node correspond to hashed bitarray of 
    # contained basic basic sets. Each node also contain attributes such as level (level), 
    # set of contained species (sp) and if the set is semi-self-maintained (ssm)
    # and an organization (is_org). 
    # The edges characterize the closures with the different basic sets, where 
    # the arrival set corresponds to the closure, the key to the 
    # basic with which the closure was performed and the attribute whether the closure 
    # is synergistic or not. 
    # The function also creates the class member (SsmStrSsmListSpArray) and (SsmStrOrgListSpArray) which 
    # corresponds to a list of the closed semi-self-maintained set and 
    # organizations respectively
    # This algorithm differs from setSynStr() by considering the basic to 
    # be conjugated which can contribute to be semi-self maintained by use of 
    # the getGBtContribBBt() function.
    def setSsmStr(self,partial_save=None,verbose=False):
        if not hasattr(self, 'GInBListBt'):
            print("The basic sets have not been initialized, please run the setGenerators() function.")
            return 
        # Initialization of the synergistic structure as a multigraph object
        G = nx.MultiDiGraph()
        # Initialization of the list of semi-self-maintained sets and organizations
        ssms=[]
        ssms_bt=[]
        org=[]
        org_bt=[]
        # step measure
        st=0
        
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.GInBListBt)):
            st+=1
            is_ssm=self.isSsmFromSp(self.BSpListBt[i])
            if is_ssm:
                is_org=self.isSmFromSp(self.BSpListBt[i])
                G.add_node(fbt(self.GInBListBt[i]),level=self.GInBListBt[i].count(),
                       sp=self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])],
                       is_ssm=is_ssm,is_org=is_org,is_basic=True,basic_id=i)
                ssms.append(self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])])
                ssms_bt.append(self.BSpListBt[i])
                if is_org:
                    org.append(self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])])
                    org_bt.append(self.BSpListBt[i])
            else:           
                G.add_node(fbt(self.GInBListBt[i]),level=self.GInBListBt[i].count(),
                       sp=self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])],
                       is_ssm=False,is_org=False,is_basic=True,basic_id=i)
        
        # Generation of the multigraph of the synergistic structure by set level (contained basics)
        for i in range(len(self.BSpListBt)):
             
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
            n_nodes=len(nodes)
            node_n=0
             # Generating closures whit connected basics sets for each set in level i
            for j in nodes:
                node_n+=1
                if verbose:
                    print("level: ",i+1, "from ",len(self.BSpListBt),", node: ",node_n," from ",n_nodes)
                # if node is semi-self-maintained (ssm), the it explore other possible 
                # combinations to search for ssm sets, if not search for basic that can 
                # contribute to be ssm
                
                if G.nodes[j]["is_ssm"]:
                     conn=self.getIndArrayFromBt(self.getGBtNotInBBt(bt(j)))
                     # connected nodes
                     conn_nodes=self.getGBtConnectedToBBt(bt(j)).search(1)
                     
                else:
                     contrib=self.getGBtContribBBt(bt(j))
                     if not contrib.any():
                         continue
                     conn=self.getIndArrayFromBt(contrib)
                     conn_nodes=contrib.search(1)

                for k in conn:
                     st+=1
                     # Closure result
                     cr_sp=self.getClosureFromSp(self.getSpBtInGBt(bt(j) | self.GInBListBt[k]),bt_type=True)
                     cr_a=fbt(self.getGBtInSpBt(cr_sp))
                     
                     # node is added if is not in structrue
                     if not (cr_a in G):
                         is_ssm=self.isSsmFromSp(cr_sp)
                         if is_ssm:
                            is_org=self.isSmFromSp(cr_sp)
                            G.add_node(cr_a,level=cr_a.count(),
                                   sp=self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)],
                                   is_ssm=is_ssm,is_org=is_org,is_basic=False)
                            ssms.append(self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)])
                            ssms_bt.append(cr_sp)
                            if is_org:
                                org.append(self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)])
                                org_bt.append(cr_sp)
                         else:           
                            G.add_node(cr_a,level=cr_a.count(),
                                   sp=self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)],
                                   is_ssm=False,is_org=False,is_basic=False)
         
                     # Adding edges corresponding to the colsure, and verifing if is a sinergy:
                     if cr_a.count() > (bt(j)|self.GInBListBt[k]).count():
                        G.add_edge(j,cr_a,key=k,syn=True,added_basic=k,conn=(k in conn_nodes))
                     else:
                        G.add_edge(j,cr_a,key=k,syn=False,added_basic=k,conn=(k in conn_nodes))
            
            if not partial_save is None:
                self.saveToPkl(partial_save) 
        
        connected_edges = [(u,v,k,e) for u,v,k,e in G.edges(keys=True,data=True) if e['conn']]
        H=nx.MultiDiGraph(connected_edges)
        nx.set_node_attributes(H,dict(G.nodes(data=True)))
        
        self.SsmStrNx=G
        self.SsmStrConnectedNx=H
        self.SsmStrSsmListSpArray=ssms
        self.SsmStrConnectedSsmListBtArray = [bt(x) for x,y in H.nodes(data=True) if y['is_ssm']]
        self.SsmStrOrgListSpArray=org
        self.SsmStrConnectedOrgListBtArray = [bt(x) for x,y in H.nodes(data=True) if y['is_org']]
        self.SsmStrSsmListBtArray=ssms_bt
        self.SsmStrOrgListBtArray=org_bt
        self.SsmStrConnectedSsmListSpArray=list(map(self.getSpBtInGBt,self.SsmStrConnectedSsmListBtArray))
        self.SsmStrConnectedOrgListSpArray=list(map(self.getSpBtInGBt,self.SsmStrConnectedOrgListBtArray))
        
        return(st)
    
    
    # Synergistic structure calculation function, requires the setGenerators() 
    # function to be executed beforehand .It returns an directed multi-graph
    # type networkx (setConnectedStr), where each node correspond to hashed bitarray of 
    # contained basic basic sets. Each node also contain attributes such as level (level), 
    # set of contained species (sp) and if the set is semi-self-maintained (ssm)
    # and an organization (is_org). 
    # The edges characterize the closures with the different basic sets, where 
    # the arrival set corresponds to the closure, the key to the 
    # basic with which the closure was performed and the attribute whether the closure 
    # is synergistic or not. 
    # The function also creates the class member (DynStrSsmListSpArray) and (DynStrOrgListSpArray) which 
    # corresponds to a list of the closed semi-self-maintained set and 
    # organizations respectively
    # This algorithm differs from setSynStr() by considering the connected basic 
    # by use of the function getGBtConnectedToBBt(), The result is a structrure 
    # where the nodes only connected closed sets.
    def setConnectedStr(self,partial_save=None,verbose=False):
        if not hasattr(self, 'GInBListBt'):
            print("The basic sets have not been initialized, please run the setGenerators() function.")
            return 
        # Initialization of the synergistic structure as a multigraph object
        G = nx.MultiDiGraph()
        # Initialization of the list of semi-self-maintained sets and organizations
        ssms=[]
        ssms_bt=[]
        org=[]
        org_bt=[]
        # number of steps
        st=0
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.GInBListBt)):
            st+=1
            is_ssm=self.isSsmFromSp(self.BSpListBt[i])
            if is_ssm:
                is_org=self.isSmFromSp(self.BSpListBt[i])
                G.add_node(fbt(self.GInBListBt[i]),level=self.GInBListBt[i].count(),
                       sp=self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])],
                       is_ssm=is_ssm,is_org=is_org,is_basic=True,basic_id=i)
                ssms.append(self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])])
                ssms_bt.append(self.BSpListBt[i])
                if is_org:
                    org.append(self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])])
                    org_bt.append(self.BSpListBt[i])
            else:           
                G.add_node(fbt(self.GInBListBt[i]),level=self.GInBListBt[i].count(),
                       sp=self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])],
                       is_ssm=False,is_org=False,is_basic=True,basic_id=i)
        
        # Generation of the multigraph of the synergistic structure by set level (contained basics)
        for i in range(len(self.BSpListBt)):
             
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
            n_nodes=len(nodes)
            node_n=0 
             # Generating closures whit connected basics sets for each set in level i
            for j in nodes:
                node_n+=1
                if verbose:
                    print("level: ",i+1, "from ",len(self.BSpListBt),", node: ",node_n," from ",n_nodes)
                # if node is semi-self-maintained (ssm), the it explore other possible 
                # combinations only whit dynamically connected basics
                # if not search for basic that can contribute to be ssm
                
                conn=self.getIndArrayFromBt(self.getGBtConnectedToBBt(bt(j)))
      
                for k in conn:
                     st+=1    
                     # Closure result
                     cr_sp=self.getClosureFromSp(self.getSpBtInGBt(bt(j) | self.GInBListBt[k]),bt_type=True)
                     cr_a=fbt(self.getGBtInSpBt(cr_sp))
                     
                     # node is added if is not in structrue
                     if not (cr_a in G):
                         is_ssm=self.isSsmFromSp(cr_sp)
                         if is_ssm:
                            is_org=self.isSmFromSp(cr_sp)
                            G.add_node(cr_a,level=cr_a.count(),
                                   sp=self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)],
                                   is_ssm=is_ssm,is_org=is_org,is_basic=False)
                            ssms.append(self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)])
                            ssms_bt.append(cr_sp)
                            if is_org:
                                org.append(self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)])
                                org_bt.append(cr_sp)
                         else:           
                            G.add_node(cr_a,level=cr_a.count(),
                                   sp=self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)],
                                   is_ssm=False,is_org=False,is_basic=False)
         
                     # Adding edges corresponding to the colsure, and verifing if is a sinergy:
                     if cr_a.count() > (bt(j)|self.GInBListBt[k]).count():
                        G.add_edge(j,cr_a,key=k,syn=True,added_basic=k)
                     else:
                        G.add_edge(j,cr_a,key=k,syn=False,added_basic=k)
            
            if not partial_save is None:
                self.saveToPkl(partial_save)    
        
        self.ConnectedStrNx=G
        self.ConnectedStrSsmListSpArray=ssms
        self.ConnectedStrOrgListSpArray=org
        self.ConnectedStrSsmListBtArray=ssms_bt
        self.ConnectedStrOrgListBtArray=org_bt
        return(st)
    
    # Synergistic structure calculation function, requires the setGenerators() 
    # function to be executed beforehand .It returns an directed multi-graph
    # type networkx (setConnectedStr), where each node correspond to hashed bitarray of 
    # contained basic basic sets. Each node also contain attributes such as level (level), 
    # set of contained species (sp) and if the set is semi-self-maintained (ssm)
    # and an organization (is_org). 
    # The edges characterize the closures with the different basic sets, where 
    # the arrival set corresponds to the closure, the key to the 
    # basic with which the closure was performed and the attribute whether the closure 
    # is synergistic or not. 
    # The function also creates the class member (DynStrSsmListSpArray) and (DynStrOrgListSpArray) which 
    # corresponds to a list of the closed semi-self-maintained set and 
    # organizations respectively
    # This algorithm differs from setSynStr() by considering the basic to 
    # be conjugated which can contribute to be semi-self maintained by use of 
    # the getGBtContribBBt() function and use of the function getGBtConnectedToBBt(), which 
    # only connect to basics if there are reactively connected. The result is a structrure 
    # where the nodes are semi-self-mantianed and only dynamically connected.
    def setSsmConnectedStr(self,partial_save=None,verbose=False):
        if not hasattr(self, 'GInBListBt'):
            print("The basic sets have not been initialized, please run the setGenerators() function.")
            return 
        # Initialization of the synergistic structure as a multigraph object
        G = nx.MultiDiGraph()
        # Initialization of the list of semi-self-maintained sets and organizations
        ssms=[]
        ssms_bt=[]
        org=[]
        org_bt=[]
        # number of steps
        st=0
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.GInBListBt)):
            st+=1
            is_ssm=self.isSsmFromSp(self.BSpListBt[i])
            if is_ssm:
                is_org=self.isSmFromSp(self.BSpListBt[i])
                G.add_node(fbt(self.GInBListBt[i]),level=self.GInBListBt[i].count(),
                       sp=self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])],
                       is_ssm=is_ssm,is_org=is_org,is_basic=True,basic_id=i)
                ssms.append(self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])])
                ssms_bt.append(self.BSpListBt[i])
                if is_org:
                    org.append(self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])])
                    org_bt.append(self.BSpListBt[i])
            else:           
                G.add_node(fbt(self.GInBListBt[i]),level=self.GInBListBt[i].count(),
                       sp=self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])],
                       is_ssm=False,is_org=False,is_basic=True,basic_id=i)
        
        # Generation of the multigraph of the synergistic structure by set level (contained basics)
        for i in range(len(self.BSpListBt)):
             
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
            n_nodes=len(nodes)
            node_n=0  
             # Generating closures whit connected basics sets for each set in level i
            for j in nodes:
                node_n+=1
                if verbose:
                    print("level: ",i+1, "from ",len(self.BSpListBt),", node: ",node_n," from ",n_nodes)
                # if node is semi-self-maintained (ssm), the it explore other possible 
                # combinations only whit dynamically connected basics
                # if not search for basic that can contribute to be ssm
                
                if G.nodes[j]["is_ssm"]:
                     conn=self.getIndArrayFromBt(self.getGBtConnectedToBBt(bt(j)))
                else:
                     contrib=self.getGBtContribBBt(bt(j))
                     if not contrib.any():
                         continue
                     conn=self.getIndArrayFromBt(contrib)
                 
                for k in conn:
                     st+=1    
                     # Closure result
                     cr_sp=self.getClosureFromSp(self.getSpBtInGBt(bt(j) | self.GInBListBt[k]),bt_type=True)
                     cr_a=fbt(self.getGBtInSpBt(cr_sp))
                     
                     # node is added if is not in structrue
                     if not (cr_a in G):
                         is_ssm=self.isSsmFromSp(cr_sp)
                         if is_ssm:
                            is_org=self.isSmFromSp(cr_sp)
                            G.add_node(cr_a,level=cr_a.count(),
                                   sp=self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)],
                                   is_ssm=is_ssm,is_org=is_org,is_basic=False)
                            ssms.append(self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)])
                            ssms_bt.append(cr_sp)
                            if is_org:
                                org.append(self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)])
                                org_bt.append(cr_sp)
                         else:           
                            G.add_node(cr_a,level=cr_a.count(),
                                   sp=self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)],
                                   is_ssm=False,is_org=False,is_basic=True)
         
                     # Adding edges corresponding to the colsure, and verifing if is a sinergy:
                     if cr_a.count() > (bt(j)|self.GInBListBt[k]).count():
                        G.add_edge(j,cr_a,key=k,syn=True,added_basic=k)
                     else:
                        G.add_edge(j,cr_a,key=k,syn=False,added_basic=k)
            
            if not partial_save is None:
                self.saveToPkl(partial_save)
                
        self.ConnectedSsmStrNx=G
        self.ConnectedSsmStrSsmListSpArray=ssms
        self.ConnectedSsmStrOrgListSpArray=org
        self.ConnectedSsmStrSsmListBtArray=ssms_bt
        self.ConnectedSsmStrOrgListBtArray=org_bt
        return(st)
    
    # Function that generates a network of the structures. It receives as 
    # input a synergic, semi-self-maintained or dynamically connected structure. 
    # The output corresponds to a pyvis object. The colors 
    # red, blue and green correspond to whether the set in question is a only 
    # reactive closed, only semi-self-maintained or an organization 
    # respectively. The shape of the set is circular if it is a basic or 
    # square if it is not. Finally, the green arrows correspond to synergies 
    # and the blue arrows to spurious union.
    def getStrDisplayPv(self,graph,notebook=False):
        G = nx.relabel_nodes(graph, lambda x: str(self.getIndArrayFromBt(bt(x))))
        nt = Network('500px', '500px',directed=True,notebook=notebook)
        # populates the nodes and edges data structures
        
        # Removing the sp array for pyvis ploting
        for i in G.nodes:
            G.nodes[i]['size']=len(G.nodes[i]['sp'])*3
            G.nodes[i]['sp']=str(G.nodes[i]['sp'])
        
        nt.from_nx(G)
        # nx.draw(G, with_labels=True, font_weight='bold')
        
        nt.toggle_physics(False)
        for i in range(len(nt.nodes)):
            # nt.nodes[i]['size']=20
            # nt.nodes[i]['title']=nt.nodes[i]['id']
            nt.nodes[i]['title']=str(nt.nodes[i]['sp'])
            if nt.nodes[i]['is_org']:
                nt.nodes[i]['color']='green'
                # nt.nodes[i]['group']=2
            elif nt.nodes[i]['is_ssm']:
                nt.nodes[i]['color']='blue'
                # nt.nodes[i]['group']=3
            else:
                nt.nodes[i]['color']='red'
                # nt.nodes[i]['group']=1
            if nt.nodes[i]['is_basic']:
                # nt.nodes[i]['title']=str(nt.nodes[i]['basic_id'])
                nt.nodes[i]['shape']='dot'
            else:
                nt.nodes[i]['shape']='square'
                
        for i in range(len(nt.edges)):
            nt.edges[i]['arrowStrikethrough']=True
            nt.edges[i]['label']=nt.edges[i]['added_basic']
            nt.edges[i]['title']=nt.edges[i]['added_basic']
            nt.edges[i]['smooth'] = True        
            if nt.edges[i]['syn']==False:
                nt.edges[i]['color']="blue"
            else:
                nt.edges[i]['color']="green"
          
        nt.directed =True  
        # nt.show('str.html')
        return(nt)
    
    # Function that returns the set directly below for a set (BtSet) of 
    # list of sets (BtList). BtSet is a bitarray object, and BtList  is a 
    # list of bitarrays.
    def getDirectlyBelowBtList(self,BtSet,BtList):
        
        # conversion for hasable bitarray object
        SortedSets=list(map(fbt,BtList))
        # Sorting elements by size (biggest to smallest)
        SortedSets=sorted(SortedSets,key= lambda x: -x.count())
        # subsetig for only contained sets 
        SortedSets=[x for x in SortedSets if subset(x, BtSet)]
        SortedSets.remove(fbt(BtSet))
        DirectlyBelowSets=[]
        while SortedSets:
            i=SortedSets[0]
            # Appendig bisgest set of the collection
            DirectlyBelowSets.append(i)
            # Removing all sets that contain currently appeded set i
            SortedSets=[x for x in SortedSets if not subset(x, i)]

            
        return DirectlyBelowSets
    
    def addAllNonReacOrgs(self,orgslist):
        '''
        

        Parameters
        ----------
        orglist : list of bitarray
            list of organization as species bitarrays.
        Returns
        -------
        Add all non reactive organization to a list of organizations.

        '''
        emptyset=bt(self.MpDf.shape[0])
        emptyset.setall(0)
        all_orgs = [com for sub in orgslist+[emptyset] for com in self.getNonReacSets(sub)]
        return(all_orgs+orgslist)
    
    
    # Generates a networkx Hasse diagram of a form a list of bitarrays 
    # Blist       
    def getHasseNxFromBtList(self,BtList,setlabel="O"):
        
        NSpSetsDict=count(list(map(lambda x: x.count(),BtList)))
        
        # conversion for hasable bitarray object
        SortedSets=list(map(fbt,BtList))
        # Sorting elements by size (smallest to biggest)
        SortedSets=sorted(SortedSets,key= lambda x: x.count())
        
        # creation of the Graph object
        Hasse=nx.Graph()
        c = 0
        node_id_count=0
        for i in SortedSets:
            size=i.count()
            
            NSpSets = NSpSetsDict[size]
            x = 0
            if NSpSets > 1:
                x = -(75*(NSpSets-1)/2) + c*(75*(NSpSets-1))
                c += 1
            else:
                c = 0
            
            y = -75*(size-1)

            Hasse.add_node(i,x=x,y=y,group=1,label=setlabel+str(node_id_count),
                           title=str(self.SpIdStrArray[self.getIndArrayFromBt(i)]),
                           fixed = json.loads('{ "x":false, "y":true}'))
            node_id_count+=1
        
        for i in SortedSets:
            
            DirectlyBelowSets=self.getDirectlyBelowBtList(i, BtList)
            if DirectlyBelowSets:
                for j in DirectlyBelowSets:
                    Hasse.add_edge(j, i,color="gray", smooth = False)
            
        return Hasse
    

    
    # Minimal generators generation function, to be started once the basic sets
    # have been generated by the setGenerators() function.  It returns a list
    # (MgenListListSpArray) of bitarray list of species that correspond to the sets so that 
    # by means of the closure they generate the basic molecule.
    def setMgen(self):
        if not hasattr(self, 'GInBListBt'):
            print("The basic sets have not been initialized, please run the setGenerators() function.")
            return 
        
        
        mgen=[]
        self.MGenStepInt=0
        # Generating a list of the support of each reaction contained in each generator
        for i in range(len(self.GRpListBt)):
            
            v=[]
            for j in self.getIndArrayFromBt(self.GRpListBt[i]):
                v.append(self.ReacListBt[j])
            
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
        self.MgenListListSpArray=mgen
    

    # Synergy generation function. it takes as input a minimum generator (sp)
    # and (pi) the index of the objective . 
    # The function verifies which combinations of generators can 
    # generate such a generator.        
    def genSyn(self,sp,pi):

        # list of generator that intersects sp 
        xp=[]     
        
        # union of generator species to see if synergy can be fulfill
        ps=bt(len(self.SpIdStrArray))
        ps.setall(0)
        # adding generator that will overlap sp
        for i in range(len(self.BSpListBt)):
               if (self.GSpListBt[i] & sp).any():
                   if(i!=pi):
                       xp.append(i)
                       ps|=self.GSpListBt[i]
        # Verifying if the generators contain the minimum sp generator 
        # and trigger the synergy.
        
        if not ((ps & sp) == sp):
            return
        
        # Bitarray for generator combinations.
        p=bt(len(xp))
        p.setall(0)
        
        # Recursive search of all synergies
        for i in range(len(xp)):
            self.recursiveGenSyn(p,i,sp,xp,pi)
            
    # Recursive synergy generation function, requires as inputs (p) the 
    # existing generators to combine, o the next level to add to the scan, 
    # xp list of indexes of the corresponding generators, (sp) minimum 
    # generator to reach and pi the index of the objective generator. 
    # Function recursively explores the possible 
    # combinations to reach a synergy. If this is reached, the branch 
    # to be explored will be cut. The synergies are stored in a list (SynReacListGBt).
    def recursiveGenSyn(self,p,o,sp,xp,pi):
        self.SynRecInt+=1
        # Species bitarray result of the sinergy
        u=bt(len(self.SpIdStrArray))
        u.setall(0)
        # Adding generator as candiadate to combinate
        p[o]=1;
        
         # Eliminating redundant generators that are already in account.
        if(p.count()>1):
            ind=self.getIndArrayFromBt(p).copy()
            for i in ind:
                p[i]=0
                u.setall(0)
                for j in self.getIndArrayFromBt(p):
                    self.MGenStepInt+=1
                    u|=self.GSpListBt[xp[j]]
                    self.SynStepsInt+=1
                p[i]=1
            
                if ((self.GSpListBt[xp[i]] & u & sp) == (self.GSpListBt[xp[i]] & sp)):
                    p[o]=0
                    return
        
        # Verfing if added generator fulfill triggering the minimal generator
        u|=self.GSpListBt[xp[o]]
        
        # If it a new synergy, it will be appned and recusrion will stop
        if ((u & sp) == sp):

            c_syn = bt(len(self.GInBListBt))
            c_syn.setall(0)
            
            for i in self.getIndArrayFromBt(p):
                self.SynStepsInt+=1
                c_syn[xp[i]]=1
            
            # if not c_syn in self.SynReacListGBt:
            self.SynReacListGBt.append(c_syn)
            op=bt(len(self.GInBListBt))
            op.setall(0)
            op[pi]=1
            self.SynProdListGBt.append(op)

                
            p[o]=0
            return
        
        # if not recursion continue
        for i in range(o+1,len(xp)):
            self.recursiveGenSyn(p,i,sp,xp,pi)
     
        p[o]=0
                
            
    # Function that generates all the synergies from the minimum 
    # generators. This is achieved through the use of the function genSyn()
    # and the recursive function recursiveGenSyn(). The output consists of list SynReacListGBt 
    # which contains all the synergies and list (SynProdListGBt) which contains 
    # all the triggered generators.   
    def setSyn(self):
        if not hasattr(self, 'MgenListListSpArray'):
            print("The minimal genetators have not been initialized, please run the setMgen(() function.")
            return         
        # List of posible synergies
        
        self.SynReacListGBt=[]
        self.SynProdListGBt=[]
        self.SynStepsInt=0
        self.SynRecInt=0
        # Generation of all synergies from all minimum generators
        for i in range(len(self.MgenListListSpArray)):
            for j in self.MgenListListSpArray[i]:
                if j.count()>0:
                    self.genSyn(j,i)
                else:
                    for k in range(len(self.MgenListListSpArray)):
                        if k!=i:
                            op=bt(len(self.GInBListBt))
                            op.setall(0)
                            op[k]=1
                            self.SynReacListGBt.append(op.copy())
                            op.setall(0)
                            op[i]=1
                            self.SynProdListGBt.append(op.copy())
                            self.SynStepsInt+=1
                    
    # Function that returns the closures that generate synergies with the 
    # input set (sp). The output corresponds to two lists in which the first 
    # corresponds to the synergy formed and the second the respective close set 
    # with which the union-closure whit (sp) set generates synergy.
    
    def getSynFromSp(self,sp_set):
        # Checks if input is or not bitarray, if it's no, it make the 
        # transformation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set

        # Only closed sets are considered, therefore sp is considered 
        # as its generated closure
        sp=self.getClosureFromSp(sp,bt_type=True)
        # generator that contain the species
        p=self.getGBtInSpBt(sp)

        # synergies that can generate synergies with p, i.e. synergies 
        # that do not contain p and but intersect p 
        syn_part=[]
        for i in self.SynReacListGBt:
            if any(p&i):
                syn_part.append(i)
        
        # list of closed set (coversion to bitarray from frozenbitarray), that 
        # not contain sp
        cl_sets=list(filter(lambda x: (not p&x==p) & (not any(p&x)),list(map(lambda x: bt(x),list(self.SynStrNx.nodes())))))

        # generation of sinergies, syn_set are the results sinergies and syn_cand
        # are the sets that generates the synergy
        syn_sets=[]
        syn_cand=[]
        for i in cl_sets:
            for j in syn_part:
                syn_p=p|i
                syn_set=self.getClosureFromSp(self.getSpBtInGBt(p|i),bt_type=True)
                if ((syn_p)&j==j) & (sp!=syn_set):
                    syn_sets.append(self.SpIdStrArray[self.getIndArrayFromBt(syn_set)])
                    syn_cand.append(self.SpIdStrArray[self.getIndArrayFromBt(self.getSpBtInGBt(i))])
                    break
                    
        return [syn_sets,syn_cand]
    
    # Function that generates the synergy interactive graph. It returns
    # pyvis object, where each synergy is labeled as p, and it's draw as a 
    # green square. The generators are colored as blue circules, where the size 
    # is proportinal to the number of contained species.
    def displaySynPv(self,notebook=False):
        
        all_part=bt(len(self.BSpListBt))
        all_part.setall(0)

        for i in range(len(self.SynReacListGBt)):
            all_part|=self.SynReacListGBt[i]
            all_part|=self.SynProdListGBt[i]

        G = nx.MultiDiGraph()

        for i in self.getIndArrayFromBt(all_part):
            sp=str(self.SpIdStrArray[self.getIndArrayFromBt(self.GSpListBt[i])])
            size=len(self.SpIdStrArray[self.getIndArrayFromBt(self.GSpListBt[i])])*3
            if self.GInBListBt[i].count()>1:
                G.add_node(i, color = "royalblue", label=str(i), size=size, title=sp, shape="dot")
            else:
                G.add_node(i, color = "lime", label=str(i), size=size, title=sp, shape="dot")
        for i in range(len(self.SynReacListGBt)):
            if self.SynReacListGBt[i].count()==1:
                G.add_node("p"+str(i), color = "crimson", label="p"+str(i), size=7, shape="square")
            else:
                G.add_node("p"+str(i), color = "teal", label="p"+str(i), size=7, shape="square")
        for i in range(len(self.SynReacListGBt)):
            if self.SynReacListGBt[i].count()==1:
                for j in self.getIndArrayFromBt(self.SynReacListGBt[i]):    
                    G.add_edge(j, "p"+str(i), color="crimson")
                for j in self.getIndArrayFromBt(self.SynProdListGBt[i]):    
                    G.add_edge("p"+str(i), j, color="crimson")
            else:        
                for j in self.getIndArrayFromBt(self.SynReacListGBt[i]):    
                    G.add_edge(j, "p"+str(i), color="teal")
                for j in self.getIndArrayFromBt(self.SynProdListGBt[i]):    
                    G.add_edge("p"+str(i), j, color="teal")
        nt = Network('500px', '500px',directed=True,notebook=notebook)
        nt.from_nx(G)
        nt.toggle_physics(False)
        # nt.show('RN.html')
        return(nt)
    
    def setAllClosedReac(self,N=None,conn_search=True,heuristic_search=True):
        
        self.AllCloseSteps=0
        
        # creation of close set lists (species and reactions)
        self.CloseReacSpBt=[]
        self.CloseReacBBt=[]
        self.CloseReacOrgSpBt=[]
        # number of basic molecules
        nb=len(self.BRpListBt)
        
        # index of generators to explore
        

        
        # case of huristic exploration
        if heuristic_search:
            pb=bt(len(self.FeasGBt.search(1)))
            
        else:
            pb=bt(nb)
            
        pb.setall(1)
        # found closed set counter
        closed_cnt=0
        
        # begining of the recursive search tree, iterating over all basic molecule
        for i in range(nb):
            #recursive search over the node
            if not self.recursiveCloseReac(pb.copy(),self.BSpListBt[i],i,
                                           closed_cnt,N,conn_search=conn_search,
                                           heuristic_search=heuristic_search):
                break
            pb[i]=0
        return closed_cnt    
        
    def recursiveCloseReac(self,pb,sp,o,closed_cnt,N,conn_search,heuristic_search):
        self.AllCloseSteps+=1
        # if explored cases exceeds N cases to explore  
        if not N is None:
            if closed_cnt>=N:
                return(0)
        
        csp=self.getClosureFromSp(sp,bt_type=True)
        ib=self.getGBtInSpBt(csp)
        
        # auxilar variable for heuristic search
        explore_sets=ib.copy()
        explore_sets.setall(0)
        
        # verifing if correspond to cases to explore
        if (not heuristic_search) and not subset(ib,pb):
            return 1
        else:
            for i in pb.search(1):
                explore_sets[self.FeasGBt.search(1)[i]]=1
            if not subset(ib,self.getClosureFromSp(self.getSpBtInGBt(explore_sets))):
                return 1
        
        # adding new closed set 
        self.CloseReacSpBt.append(csp)
        self.CloseReacBBt.append(ib)
        if self.isSsmFromSp(csp):
            if self.isSsmFromSp(csp):
                self.CloseReacOrgSpBt.append(csp)
        closed_cnt+=1
        o+=1
        
        # if conn_search:
        #     pb&=self.getGBtConnectedToBBt(ib,inculde_itself=True)
        # # if heuristic_search:
        #     pb&=self.getGBtFeasbleConnectedtoBBt(ib,True)
        iter_range=np.array(pb.search(1))
        iter_range=iter_range[iter_range>=o]
        print_range=iter_range.copy()
        explore_sets=ib.copy()
        
        for i in iter_range:
            
            if heuristic_search:
                ind = self.FeasGBt.search(1)[i]
            else:
                ind=i
            print("exploting cases",print_range)
            if ib[i]==1:
                print_range=print_range[print_range!=i]
                continue
            
            add_gen=ib.copy()
            add_gen.setall(0)
            add_gen[ind]=1
            print("current added generator",self.SpIdStrArray[self.getSpBtInGBt(add_gen).search(1)])
            print("current closed set",self.SpIdStrArray[self.getSpBtInGBt(ib).search(1)])    
            print("connectivity condition",subset(add_gen,self.getGBtTransConnectedToBBt(ib)))
            print("heuristic condition",subset(add_gen,self.FeasGBt))
            
            if heuristic_search:
                explore_sets.setall(0)
                for i in pb.search(1):
                    try:
                        explore_sets[self.FeasGBt.search(1)[i]]=1
                    except:
                        print("i",i)
                        
                condition=(explore_sets&self.GInBListBt[ind]).any() # and not (ib&self.GInBListBt[ind]).any()
            else:
                condition=subset(self.GInBListBt[i],pb)
               
            if condition:
                if (not heuristic_search) and (not conn_search):
                    
                    new_set=csp|self.BSpListBt[ind]
                    if not self.recursiveCloseReac(pb.copy(), new_set, i, closed_cnt, 
                                                    N,conn_search=conn_search,heuristic_search=heuristic_search):
                        return 0
                
                elif heuristic_search and (not conn_search) and subset(add_gen,self.FeasGBt):
                    print("entre")
                    new_set=csp|self.BSpListBt[ind]
                    if not self.recursiveCloseReac(pb.copy(), new_set, i, closed_cnt, 
                                                    N,conn_search=conn_search,heuristic_search=heuristic_search):
                        return 0

                elif (not heuristic_search) and conn_search and subset(add_gen,self.getGBtTransConnectedToBBt(ib)):
                    # print("current added generator",RN.SpIdStrArray[add_gen.search(1)])
                    # print("current closed set",RN.SpIdStrArray[ib.search(1)])
                    new_set=csp|self.BSpListBt[ind]
                    if not self.recursiveCloseReac(pb.copy(), new_set, i, closed_cnt, 
                                                   N,conn_search=conn_search,heuristic_search=heuristic_search):
                        return 0
                    
                elif heuristic_search and conn_search and subset(add_gen,self.getGBtTransConnectedToBBt(ib)) and\
                    subset(add_gen,self.FeasGBt):
                    
                    new_set=csp|self.BSpListBt[ind]
                    if not self.recursiveCloseReac(pb.copy(), new_set, i, closed_cnt, 
                                                    N,conn_search=conn_search,heuristic_search=heuristic_search):
                        return 0
            pb[i]=0
            
        return 1
        
    def isFeasibleGBt(self,gen_ind):
                
        # stoichiometric matrix of rn with prepended destruction reactions
        S=self.MpDf-self.MrDf
        M=100*self.MpDf.shape[0]*self.MpDf.shape[1]
        # Name and type of problem
        prob = LpProblem("LP_Org_test", LpMinimize)
        x={}
        for i in range(self.MpDf.shape[1]): 

            if i in self.GRpListBt[gen_ind].search(1): 
                x_i = LpVariable(f"process_{i}", lowBound=1, cat='Continuous')
            else:
                x_i = LpVariable(f"process_{i}", lowBound=0, cat='Continuous')
            
            
            
            x[i]=x_i
            
        y = LpVariable.dicts("upbond", S.iloc[0].index, lowBound=1, cat='integer')
        # We define the objective function
        
        objective = lpSum([x[i] for i in x])+M*lpSum([y[i] for i in y])
        prob += objective
       
        # prob += lpSum([x[i] for i in x])
        for j in S.index:
            prob +=lpSum([S.loc[j][i] * x[i] for i in S.loc[j].index]) >= 0
            prob += x[i] <= y[i]*M
        # Resolver el problema
        prob.solve(solver=GLPK(msg=False))

        if LpStatus[prob.status]=='Optimal':
            v_result=np.zeros(self.MpDf.shape[1])
            
            for v in prob.variables():
                
                if ("process" in v.name):

                    ind=int(re.search(r'\d+', v.name).group())
                    v_result[ind]=v.varValue
            
            reac_ind=np.where(v_result>0)[0]
            reac_ind=list((set(reac_ind)-set(self.GRpListBt[gen_ind].search(1))))

            reac_bt=bt(self.MpDf.shape[1])
            reac_bt.setall(0)
            
            if len(reac_ind)>0:
                for i in reac_ind:
                    reac_bt[i]=1
                
            required_G=bt(len(self.GRpListBt))
            required_G.setall(0)
            for ind, i in enumerate(self.GRpListBt):
                if (i&reac_bt).any():
                    required_G[ind]=1
                    
            return True, v_result, required_G
        else:
            
            required_G=bt(len(self.GRpListBt))
            required_G.setall(0)
            
            return False, np.zeros(self.MpDf.shape[1]), required_G
        
        
    def setFeasG(self):
        
        feasible_generators=[]
        self.ReqGBtforFeas=[]
        
        for i in range(len(self.GRpListBt)):
            feasible_generators.append(self.isFeasibleGBt(i))
            self.ReqGBtforFeas.append(feasible_generators[-1][2])
            
        G_candidates=bt(len(feasible_generators))
        G_candidates.setall(1)
 
        NullGbt=G_candidates.copy()
        NullGbt.setall(0)
        
        remcand=True

        while(remcand):
            remcand=False
            for i in G_candidates.search(1):
                if feasible_generators[i][0]==False:
                    G_candidates[i]=0
                    remcand=True
                if not subset(feasible_generators[i][2],G_candidates):
                    G_candidates[i]=0
                    remcand=True
                    
        self.FeasGBt=G_candidates
        
        self.ReqGIndforFeas=G_candidates.search(1)
        
        G_to_remove=G_candidates.search(0)
        G_to_remove.sort(reverse=True)
        
        for i in G_to_remove:
            self.ReqGBtforFeas.pop(i)
        
        useful_reac=bt(self.MpDf.shape[1])
        useful_reac.setall(0)
        for i in self.FeasGBt.search(1):
            useful_reac|=self.GRpListBt[i]
        
        self.FeasRpBt=useful_reac
        
        G_to_remove=[]
        self.FeasHelp2SustSp=[]
        k=0
        for i in self.FeasGBt.search(1):
            S_R_dict= self.getHelp2SustVer(self.BSpListBt[i])
            # S_R_bt = reduce(lambda x,y : x[1]|y[1], self.getHelp2SustVer(self.BSpListBt[i]))
            S_R_bt=bt(self.MpDf.shape[0])
            S_R_bt.setall(0)
            
            for j in S_R_dict.keys():
                    S_R_bt|=S_R_dict[j]
            
            if not (S_R_bt&self.GSpListBt[i]).any():
                G_to_remove.append(k)
                self.FeasGBt[i]=0
                
            k+=1
            
        G_to_remove.sort(reverse=True)
        for i in G_to_remove:
            self.ReqGBtforFeas.pop(i)
        
        self.ReqGIndforFeas=self.FeasGBt.search(1)
        
        useful_reac=bt(self.MpDf.shape[1])
        useful_reac.setall(0)
        for i in self.FeasGBt.search(1):
            useful_reac|=self.GRpListBt[i]
        
        self.FeasRpBt=useful_reac
        
    def getHelp2SustVer(self,sp_set,r_set=None,useful_reac=False):
        
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
            n_reac = r.search(1)
        
        n_reac=np.array(list(set(n_reac)&set(self.getTriggerableRpBtFromSp(sp_set).search(1))))
        if useful_reac:
            n_reac=np.array(list(set(n_reac)&set(self.FeasRpBt.search(1))))
        # creating acutocatalityc set sustainable sets and the help to sustain 
        # set
        S_R_dict=dict()
        
        for i in n_reac:
            S_R_dict[i]=bt(self.MpDf.shape[0])
            S_R_dict[i].setall(0)
        # assigning the variables to the sets described above
        for i in n_reac:
        
            sp_i=self.ReacListBt[i]|self.ProdListBt[i]
            
            # estimating which of the reacting species are autocatalytic and sustained
            auto_cata=(self.MpDf.iloc[sp_i.search(1),i]-self.MrDf.iloc[sp_i.search(1),i])>0
            auto_cata=auto_cata[auto_cata]
            auto_cata=np.where(self.SpIdStrArray == np.array(auto_cata.index)[:,None])[1]
            sustain=(self.MpDf.iloc[sp_i.search(1),i]-self.MrDf.iloc[sp_i.search(1),i])>=0
            sustain=sustain[sustain]
            sustain=np.where(self.SpIdStrArray == np.array(sustain.index)[:,None])[1]
            
            a_i=bt(self.MpDf.shape[0])
            s_i=a_i.copy()
            a_i.setall(0)
            s_i.setall(0)
            
            # assinging variables of the current reaction
            for k in sustain:
                s_i[k]=1
                if k in auto_cata:
                    a_i[k]=1
            
            if self.ReacListBt[i].any():
                S_R_dict[i]|=s_i&self.ReacListBt[i]
            else:
                S_R_dict[i]|=s_i
            
            n_reac_j=n_reac.copy()
            n_reac_j=np.delete(n_reac_j,np.where(n_reac==i)[0])
            
            for j in n_reac_j:
                
                sp_j=self.ReacListBt[j]|self.ProdListBt[j]
                sustain=(self.MpDf.iloc[sp_j.search(1),j]-self.MrDf.iloc[sp_j.search(1),j])>=0
                sustain=sustain[sustain]
                sustain=np.where(self.SpIdStrArray == np.array(sustain.index)[:,None])[1]
                s_j=s_i.copy()
                s_j.setall(0)
                for k in sustain:
                    s_j[k]=1

                if s_j.count(1)>0 and (self.ReacListBt[j]&a_i).count(1)>0: 
                    
                    S_R_dict[i]|=s_i&self.ReacListBt[j]

        return(S_R_dict)
    
                
    

    # For a given bitarray of contained basic set, return the dynamically 
    # connected basics sets to it
    def getGBtFeasbleConnectedtoBBt(self,s,inculde_itself=False):
        c=bt(len(self.BSpListBt))
        c.setall(0)
                
        for i in self.FeasGBt.search(1):
            j=np.where(np.array(self.FeasGBt.search(1))==i)[0][0]
            if (s&self.ReqGBtforFeas[j]).any():
                c[i]=1
    
        if not inculde_itself:
            c= c & ~s
        else:
            c= c | s
        return c
        
    def getBoolExploreFactofBBt(self,s,i):

        try:
            j=np.where(np.array(self.FeasGBt.search(1))==i)[0][0]
            if not subset(self.ReqGBtforFeas[j],s):
                return False
        except:
            return False
        
        return True
        
        
            
        
                