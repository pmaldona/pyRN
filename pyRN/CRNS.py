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

# Class that calculates the synergistic and organizational structure 
# of the closed reactants of a reaction network. 
class CRNS(RNIRG):
    
    # Function that generates the basic sets of a reaction network. When 
    # executed, it generates the following members of the class, which 
    # correspond to a list of bitarrays, where each element of the list 
    # corresponds to a basic set or partition as the case may be:
    
    # sp_b: species contained in each basic set
    # r_b: reactions contained in each basic set 
    # sp_p: species contained in each partition
    # r_p: reactions contained in each partition 
    # rsp_b: species of reactants contained in each basic 
    # psp_b: species of products contained in each basic
    # sp_sp_b: positive stoichiometric species contained in the closure of each partition
    # sn_sp_b: negative stoichiometric species contained in the closure of each partition
    # p_b: partitions contained in each basic set
    # dyn_conn: basic sets that are connected (whose has a non empty intersection 
    # between the support and the products o viceversa).
            
    # The variable x_r_p corresponds to a vector of integers, of the length 
    # of the number of reactions, where the integer indicates to which 
    # partition the reaction belongs.
    def gen_basics(self):
        # creating all closed set for closure of reactant part
        c_reac=[]
          
        for i in self.reac:
            c_reac.append(self.closure(i,bt_type=True))
        
        # number of equivalence clases
        xeqc=0
        # class equivalence vector
        x_r_p=-np.ones(self.mp.shape[1])
        x_r_p=x_r_p.tolist()
        # list of basic sets
        sp_b=[]
        
        # number of steps
        st=0
        # creating equivalence classes for each reaction that generate the same
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
        
        
        # assigning reaction contained each partition
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
        r_b=[] # reactions supported by each basic (equivalence class)
        rsp_b=[] # reactants species contained each basic (equivalence class)
        psp_b=[] # product species contained each basic (equivalence class)
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
       
        # Creating connecting basics and dynamically connected basics sets
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
    
    # Function that returns the number of basic sets in which each species of 
    # the vector sp appears. The input vector sp can be a list of 
    # strings or a bitarray.
    def basic_sp_presence(self,sp_set):
        # Checks if input is or not bitarray, if it's no, it make the 
        # transformation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
        
        sp_b_presc=np.zeros(len(sp))
        sp_p_presc=np.zeros(len(sp))
        for i in sp.itersearch(1):
            for j in range(len(self.sp_b)):
                if self.sp_b[j][i]==1:
                    sp_b_presc[i]+=1
            for j in range(len(self.sp_p)):
                if self.sp_p[j][i]==1:
                    sp_p_presc[i]+=1
        sp_b_presc=sp_b_presc[self.bt_ind(sp)]
        sp_p_presc=sp_p_presc[self.bt_ind(sp)]
        
        sp_b_presc=pd.DataFrame(sp_b_presc)
        sp_b_presc=sp_b_presc.transpose()
        sp_b_presc.columns=self.sp[self.bt_ind(sp)]
        
        sp_p_presc=pd.DataFrame(sp_p_presc)
        sp_p_presc=sp_p_presc.transpose()
        sp_p_presc.columns=self.sp[self.bt_ind(sp)]
        
        return([sp_b_presc,sp_p_presc])
    
    
    # Function that returns the number of basic sets in which each reaction of 
    # the vector r appears. The input vector v can be a list of 
    # strings or a bitarray.
    def basic_r_presence(self,r_set):
        # Checks if input is or not bitarray, if it's no, it make the 
        # transformation
        if not (isinstance(r_set,bt)):
            v=bt(self.mp.shape[1])
            v.setall(0)
            
            for i in r_set:
                v[i]=1
        else:
            v=r_set

        
        r_b_presc=np.zeros(len(v))
        
        for i in v.itersearch(1):
            for j in range(len(self.r_b)):
                if self.r_b[j][i]==1:
                    r_b_presc[i]+=1

        r_b_presc=r_b_presc[self.bt_ind(v)]
        
        r_b_presc=pd.DataFrame(r_b_presc)
        r_b_presc=r_b_presc.transpose()
        
        names=[]
        
        for i in self.mp.columns:
            names.append("r"+str(i))
        
        r_b_presc.columns=np.array(names)[self.bt_ind(v)]
          
        return(r_b_presc)
    
    # Function that returns the number of basic sets in which each species of 
    # the vector sp appears. The input vector sp can be a list of 
    # strings or a bitarray.
    def plot_basic_sp_presence(self,sp_set):
        histo=self.basic_sp_presence(sp_set)
        
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
        
        histo=self.basic_sp_presence(sp_set)
        
        names=histo[1].columns.to_list()
        values=histo[1].iloc[0].tolist()
    
        fig = plt.figure(figsize = (10, 5))
 
        # creating the bar plot
        plt.bar(names, values, color ='maroon',
                width = 0.4)
         
        plt.xlabel("Species")
        plt.ylabel("Number of partitions")
        plt.title("Number of partitions sets that contain each species")
        plt.show()
    
    
    # Function that plots the number of basic sets in which each reaction of 
    # the vector r appears. The input vector v can be a list of 
    # strings or a bitarray.
    def plot_basic_r_presence(self,v):
        
        histo=self.basic_r_presence(v)
        
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
    
    
    # Function of species that returns the partition that contains the species
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
    
    
    # For a given bitarray of contained basic set, return the connected basics 
    # (basics that are not contained)
    def conn_b(self,s):
        
        c=bt(len(self.sp_b))
        c.setall(0)
        
        for i in self.bt_ind(s):
            c|=self.conn[i]
            
        c= c & ~s
        return c
        
    
    # For a given bitarray of contained basic set, return the dynamically 
    # connected basics sets to it
    def dyn_conn_b(self,s):
        c=bt(len(self.sp_b))
        c.setall(0)
        
        for i in self.bt_ind(s):
            c|=self.dyn_conn[i]
            
        c= c & ~s
        return c
    
    
    # For a given bitarray of contained basic set, return the basic sets that 
    # can contribute for the current set to be semi-self-maintained 
    def contrib_b(self,s):
        
        # negative sotichiometric species
        
        n=bt(self.mp.shape[0])
        n.setall(0)
        p=n.copy()
        pp=n.copy()
        
        for i in self.bt_ind(s):
            p|=self.sp_sp_b[i]
            n|=self.sn_sp_b[i]
                       
        # p correspond to bitarray whit positive stoichiometry
        n= n & ~p
        
        c=s.copy()
        if not n.any():
           c.setall(0)
           return c
       
        # possible sets that can contribute to be semi-self-maintained
        c = ~c
        
        # eliminating basics form c that no will contribute
        for i in self.bt_ind(c):
            if not (self.sp_sp_b[i] & n).any():
                c[i]=0
            else:
                pp|=self.sp_sp_b[i]
        
        # pp are possible produces species if they do not fulfill consumed 
        # species, then there is no contribution to be a semi-self-maintained

        if (n & ~pp).any():
            c.setall(0)
            return c
        else:
            return c
        
    
    # Synergistic structure calculation function, requires the gen_basics() 
    # function to be executed beforehand .It returns an directed multi-graph
    # type networkx (syn_str), where each node correspond to hashed bitarray of 
    # contained basic basic sets. Each node also contain attributes such as level (level), 
    # set of contained species (sp) and if the set is semi-self-maintained (ssm)
    # and an organization (is_org). 
    # The edges characterize the closures with the different basic sets, where 
    # the arrival set corresponds to the closure, the key to the 
    # basic with which the closure was performed and the attribute whether the closure 
    # is synergistic or not. 
    # The function also creates the class member (syn_ssm) and (syn_org) which 
    # corresponds to a list of the closed semi-self-maintained set and 
    # organizations respectively
    def gen_syn_str(self,org_cal=False):
        if not hasattr(self, 'p_b'):
            print("The basic sets have not been initialized, please run the gen_basics() function.")
            return 
        # Initialization of the synergistic structure as a multigraph object
        G = nx.MultiDiGraph()
        # Initialization of the list of semi-self-maintained sets and organizations
        ssms=[]
        org=[]
        # step measure
        st=0
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.p_b)):
            st+=1
            is_ssm=self.is_ssm(self.sp_b[i])
            if is_ssm:
                is_org=self.is_ssm(self.sp_b[i])
                G.add_node(fbt(self.p_b[i]),level=self.p_b[i].count(),
                       sp=self.sp[self.bt_ind(self.sp_b[i])],
                       is_ssm=is_ssm,is_org=is_org,is_basic=True,basic_id=i)
                ssms.append(self.sp[self.bt_ind(self.sp_b[i])])
                if is_org:
                    org.append(self.sp[self.bt_ind(self.sp_b[i])])
            else:           
                G.add_node(fbt(self.p_b[i]),level=self.p_b[i].count(),
                       sp=self.sp[self.bt_ind(self.sp_b[i])],
                       is_ssm=False,is_org=False,is_basic=True,basic_id=i)

                
        # Generation of the multigraph of the synergistic structure by set level (contained basics)
        for i in range(len(self.sp_b)):
             
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
             
             # Generating closures whit connected basics sets for each set in level i
            for j in nodes:
                 for k in self.bt_ind(self.conn_b(bt(j))):
                     st+=1    
                     # Closure result
                     cr_sp=self.closure(self.p2sp(bt(j) | self.p_b[k]),bt_type=True)
                     cr_p=fbt(self.sp2p(cr_sp))
                     
                     # node is added if is not in structrue
                     if not (cr_p in G):
                         is_ssm=self.is_ssm(cr_sp)
                         if is_ssm:
                             is_org=self.is_ssm(cr_sp)
                             G.add_node(cr_p,level=cr_p.count(),
                                    sp=self.sp[self.bt_ind(cr_sp)],
                                    is_ssm=is_ssm,is_org=is_org,is_basic=False)
                             ssms.append(self.sp[self.bt_ind(cr_sp)])
                             if is_org:
                                 org.append(self.sp[self.bt_ind(cr_sp)])
                         else:           
                             G.add_node(cr_p,level=cr_p.count(),
                                    sp=self.sp[self.bt_ind(cr_sp)],
                                    is_ssm=False,is_org=False,is_basic=False)
                         
                     # Adding edges corresponding to the colsure, and verifing if is a sinergy:
                     if cr_p.count() > (bt(j)|self.p_b[k]).count():
                        G.add_edge(j,cr_p,key=fbt(self.p_b[k]),syn=True,added_basic=k)
                     else:
                        G.add_edge(j,cr_p,key=fbt(self.p_b[k]),syn=False,added_basic=k)
                    # if cr_p.count() > (bt(j)|self.p_b[k]).count():
                    #    G.add_edge(j,cr_p,key=fbt(self.p_b[k]),syn=True,added_basic=k)
                    # else:
                    #    G.add_edge(j,cr_p,key=fbt(self.p_b[k]),syn=False,added_basic=k)
                
        self.syn_str=G
        self.syn_ssms=ssms
        self.syn_org=org
        return(st)
                    

    # Synergistic structure calculation function, requires the gen_basics() 
    # function to be executed beforehand .It returns an directed multi-graph
    # type networkx (ssm_str), where each node correspond to hashed bitarray of 
    # contained basic basic sets. Each node also contain attributes such as level (level), 
    # set of contained species (sp) and if the set is semi-self-maintained (ssm)
    # and an organization (is_org). 
    # The edges characterize the closures with the different basic sets, where 
    # the arrival set corresponds to the closure, the key to the 
    # basic with which the closure was performed and the attribute whether the closure 
    # is synergistic or not. 
    # The function also creates the class member (ssm_ssms) and (ssm_org) which 
    # corresponds to a list of the closed semi-self-maintained set and 
    # organizations respectively
    # This algorithm differs from ge_syn_str() by considering the basic to 
    # be conjugated which can contribute to be semi-self maintained by use of 
    # the contrib_b() function.
    def gen_ssm_str(self):
        if not hasattr(self, 'p_b'):
            print("The basic sets have not been initialized, please run the gen_basics() function.")
            return 
        # Initialization of the synergistic structure as a multigraph object
        G = nx.MultiDiGraph()
        # Initialization of the list of semi-self-maintained sets and organizations
        ssms=[]
        org=[]
        # step measure
        st=0
        
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.p_b)):
            st+=1
            is_ssm=self.is_ssm(self.sp_b[i])
            if is_ssm:
                is_org=self.is_ssm(self.sp_b[i])
                G.add_node(fbt(self.p_b[i]),level=self.p_b[i].count(),
                       sp=self.sp[self.bt_ind(self.sp_b[i])],
                       is_ssm=is_ssm,is_org=is_org)
                ssms.append(self.sp[self.bt_ind(self.sp_b[i])])
                if is_org:
                    org.append(self.sp[self.bt_ind(self.sp_b[i])])
            else:           
                G.add_node(fbt(self.p_b[i]),level=self.p_b[i].count(),
                       sp=self.sp[self.bt_ind(self.sp_b[i])],
                       is_ssm=False,is_org=False)
        
        # Generation of the multigraph of the synergistic structure by set level (contained basics)
        for i in range(len(self.sp_b)):
             
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
             
             # Generating closures whit connected basics sets for each set in level i
            for j in nodes:
                # if node is semi-self-maintained (ssm), the it explore other possible 
                # combinations to search for ssm sets, if not search for basic that can 
                # contribute to be ssm
                
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
                     cr_sp=self.closure(self.p2sp(bt(j) | self.p_b[k]),bt_type=True)
                     cr_p=fbt(self.sp2p(cr_sp))
                     
                     # node is added if is not in structrue
                     if not (cr_p in G):
                         is_ssm=self.is_ssm(cr_sp)
                         if is_ssm:
                            is_org=self.is_ssm(cr_sp)
                            G.add_node(cr_p,level=cr_p.count(),
                                   sp=self.sp[self.bt_ind(cr_sp)],
                                   is_ssm=is_ssm,is_org=is_org)
                            ssms.append(self.sp[self.bt_ind(cr_sp)])
                            if is_org:
                                org.append(self.sp[self.bt_ind(cr_sp)])
                         else:           
                            G.add_node(cr_p,level=cr_p.count(),
                                   sp=self.sp[self.bt_ind(cr_sp)],
                                   is_ssm=False,is_org=False)
         
                     # Adding edges corresponding to the colsure, and verifing if is a sinergy:
                     if cr_p.count() > (bt(j)|self.p_b[k]).count():
                        G.add_edge(j,cr_p,key=k,syn=True)
                     else:
                        G.add_edge(j,cr_p,key=k,syn=False)
                
        self.ssm_str=G
        self.ssm_ssms=ssms
        self.ssm_org=org
        return(st)
    
    
    # Synergistic structure calculation function, requires the gen_basics() 
    # function to be executed beforehand .It returns an directed multi-graph
    # type networkx (ssm_str), where each node correspond to hashed bitarray of 
    # contained basic basic sets. Each node also contain attributes such as level (level), 
    # set of contained species (sp) and if the set is semi-self-maintained (ssm)
    # and an organization (is_org). 
    # The edges characterize the closures with the different basic sets, where 
    # the arrival set corresponds to the closure, the key to the 
    # basic with which the closure was performed and the attribute whether the closure 
    # is synergistic or not. 
    # The function also creates the class member (ssm_ssms) and (ssm_org) which 
    # corresponds to a list of the closed semi-self-maintained set and 
    # organizations respectively
    # This algorithm differs from ge_syn_str() by considering the basic to 
    # be conjugated which can contribute to be semi-self maintained by use of 
    # the contrib_b() function and use of the function d_connect(), which 
    # only connect to basics if there are reactively connected. The result is a structrure 
    # where the nodes are semi-self-mantianed and only dynamically connected.
    def gen_dyn_str(self):
        if not hasattr(self, 'p_b'):
            print("The basic sets have not been initialized, please run the gen_basics() function.")
            return 
        # Initialization of the synergistic structure as a multigraph object
        G = nx.MultiDiGraph()
        # Initialization of the list of semi-self-maintained sets and organizations
        ssms=[]
        org=[]
        # number of steps
        st=0
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.p_b)):
            st+=1
            is_ssm=self.is_ssm(self.sp_b[i])
            if is_ssm:
                is_org=self.is_ssm(self.sp_b[i])
                G.add_node(fbt(self.p_b[i]),level=self.p_b[i].count(),
                       sp=self.sp[self.bt_ind(self.sp_b[i])],
                       is_ssm=is_ssm,is_org=is_org)
                ssms.append(self.sp[self.bt_ind(self.sp_b[i])])
                if is_org:
                    org.append(self.sp[self.bt_ind(self.sp_b[i])])
            else:           
                G.add_node(fbt(self.p_b[i]),level=self.p_b[i].count(),
                       sp=self.sp[self.bt_ind(self.sp_b[i])],
                       is_ssm=False,is_org=False)
        
        # Generation of the multigraph of the synergistic structure by set level (contained basics)
        for i in range(len(self.sp_b)):
             
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
             
             # Generating closures whit connected basics sets for each set in level i
            for j in nodes:
                
                # if node is semi-self-maintained (ssm), the it explore other possible 
                # combinations only whit dynamically connected basics
                # if not search for basic that can contribute to be ssm
                
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
                     cr_sp=self.closure(self.p2sp(bt(j) | self.p_b[k]),bt_type=True)
                     cr_p=fbt(self.sp2p(cr_sp))
                     
                     # node is added if is not in structrue
                     if not (cr_p in G):
                         is_ssm=self.is_ssm(cr_sp)
                         if is_ssm:
                            is_org=self.is_ssm(cr_sp)
                            G.add_node(cr_p,level=cr_p.count(),
                                   sp=self.sp[self.bt_ind(cr_sp)],
                                   is_ssm=is_ssm,is_org=is_org)
                            ssms.append(self.sp[self.bt_ind(cr_sp)])
                            if is_org:
                                org.append(self.sp[self.bt_ind(cr_sp)])
                         else:           
                            G.add_node(cr_p,level=cr_p.count(),
                                   sp=self.sp[self.bt_ind(cr_sp)],
                                   is_ssm=False,is_org=False)
         
                     # Adding edges corresponding to the colsure, and verifing if is a sinergy:
                     if cr_p.count() > (bt(j)|self.p_b[k]).count():
                        G.add_edge(j,cr_p,key=k,syn=True)
                     else:
                        G.add_edge(j,cr_p,key=k,syn=False)
                
        self.dyn_ssm_str=G
        self.dyn_ssms=ssms
        self.dyn_org=org
        return(st)
    
    # Function that generates a network of the structures. It receives as 
    # input a synergic, semi-self-maintained or dynamically connected structure. 
    # The output corresponds to a pyvis object. The colors 
    # red, blue and green correspond to whether the set in question is a only 
    # reactive closed, only semi-self-maintained or an organization 
    # respectively. The shape of the set is circular if it is a basic or 
    # square if it is not. Finally, the green arrows correspond to synergies 
    # and the blue arrows to spurious union.
    def display_str(self,graph):
        G = nx.relabel_nodes(graph, lambda x: str(self.bt_ind(bt(x))))
        nt = Network('500px', '500px',directed=True)
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
    

    # Proto-synergy generation function. it takes as input a minimum generator (sp)
    # and (pi) the index of the objective partition. 
    # The function verifies which combinations of partitions can 
    # generate such a generator.        
    def syn_gen(self,sp,pi):

        # list of partition that intersects sp 
        xp=[]     
        
        # union of partition species to see if synergy can be fulfill
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
        
            
    # Recursive proto synergy generation function, requires as inputs (p) the 
    # existing partitions to combine, o the next level to add to the scan, 
    # xp list of indexes of the corresponding partitions, (sp) minimum 
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
    # which contains all the proto synergies and list (syn_p) which contains 
    # all the triggered partitions.   
    def all_syn(self):
        if not hasattr(self, 'mgen'):
            print("The minimal genetators have not been initialized, please run the gen_mgen() function.")
            return         
        # List of posible synergies
        
        self.syn=[]
        self.syn_p=[]
        
        # Generation of all proto synergies from all minimum generators
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
      
                    
    # Function that returns the closures that generate synergies with the 
    # input set (sp). The output corresponds to two lists in which the first 
    # corresponds to the synergy formed and the second the respective close set 
    # with which the union-closure whit (sp) set generates synergy.
    def syn_sets(self,sp_set):
        # Checks if input is or not bitarray, if it's no, it make the 
        # transformation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.index.values:
                    ind=self.mp.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set

        # Only closed sets are considered, therefore sp is considered 
        # as its generated closure
        sp=self.closure(sp,bt_type=True)
        # partition that contain the species
        p=self.sp2p(sp)

        # Proto-synergies that can generate synergies with p, i.e. synergies 
        # that do not contain p and but intersect p 
        syn_part=[]
        for i in self.syn:
            if any(p&i):
                syn_part.append(i)
        
        # list of closed set (coversion to bitarray from frozenbitarray), that 
        # not contain sp
        cl_sets=list(filter(lambda x: (not p&x==p) & (not any(p&x)),list(map(lambda x: bt(x),list(self.syn_str.nodes())))))

        # generation of sinergies, syn_set are the results sinergies and syn_cand
        # are the sets that generates the synergy
        syn_sets=[]
        syn_cand=[]
        for i in cl_sets:
            for j in syn_part:
                syn_p=p|i
                syn_set=self.closure(self.p2sp(p|i),bt_type=True)
                if ((syn_p)&j==j) & (sp!=syn_set):
                    syn_sets.append(self.sp[self.bt_ind(syn_set)])
                    syn_cand.append(self.sp[self.bt_ind(self.p2sp(i))])
                    break
                    
        return [syn_sets,syn_cand]
    
    # Function that generates the proto-synergiy interactive graph. It returns
    # oyvis objecto. whre each proto-synergie is label as p, and in draw as a 
    # green squere. The partitions are colored as blue circules, where the size 
    # is proportinal to the number of contained species.
    def display_pr_syn(self):
        
        all_part=bt(len(self.sp_b))
        all_part.setall(0)

        for i in range(len(self.syn)):
            all_part|=self.syn[i]
            all_part|=self.syn_p[i]

        G = nx.MultiDiGraph()

        for i in self.bt_ind(all_part):
            sp=str(self.sp[self.bt_ind(self.sp_p[i])])
            size=len(self.sp[self.bt_ind(self.sp_p[i])])*3
            G.add_node(i, color = "blue", label=str(i), size=size, title=sp, shape="dot")
            
        for i in range(len(self.syn)):
            G.add_node("p"+str(i), color = "green", label="p"+str(i), size=7, shape="square")
            
        for i in range(len(self.syn)):
            for j in self.bt_ind(self.syn[i]):    
                G.add_edge(j, "p"+str(i), color="gray")
            for j in self.bt_ind(self.syn_p[i]):    
                G.add_edge("p"+str(i), j, color="gray")
            
        nt = Network('500px', '500px',directed=True)
        nt.from_nx(G)
        nt.toggle_physics(False)
        # nt.show('RN.html')
        return(nt)