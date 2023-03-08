#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 20:12:37 2022

@author: pmaldona

Closed Reaction Network Strucutre Library
"""
from .CRNS import CRNS
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
from joblib import Parallel, delayed

# Class that calculates the synergistic and organizational structure 
# of the closed reactants of a reaction network. 
class CRNSMP(CRNS):
    
    def addGElements(self,i):
        '''
        

        Parameters
        ----------
        i : int
            iterator for the genertors.

        Returns
        -------
        returns elements of generators properties. Function implementend 
        for multicore parallelization.

        '''
        b_sp=bt(self.MpDf.shape[0])
        b_sp.setall(0)
        b_r=bt(self.MpDf.shape[1])
        b_r.setall(0)
        
        GInBListBt=self.getGBtInSpBt(self.BSpListBt[i])
        BRpListBt=b_r.copy()
        BReacSpListBt=b_sp.copy()
        BProdSpListBt=b_sp.copy()
        BStoichioPositiveSpListBt=b_sp.copy()
        BStoichioNegativeSpListBt=b_sp.copy()
       
        for j in GInBListBt.search(1):
            BRpListBt|=self.GRpListBt[j]
            for k in self.GRpListBt[j].search(1): 
                BReacSpListBt|=self.ReacListBt[k]
                BProdSpListBt|=self.ProdListBt[k]
                BStoichioPositiveSpListBt|=bt(self.MrDf.iloc[:,k] < self.MpDf.iloc[:,k])
                BStoichioNegativeSpListBt|=bt(self.MrDf.iloc[:,k] > self.MpDf.iloc[:,k])
        return([GInBListBt,BRpListBt,BReacSpListBt,BProdSpListBt,BStoichioPositiveSpListBt,BStoichioNegativeSpListBt])
    
    def connGGen(self,i):

        b_c=bt(len(self.BSpListBt))
        b_c.setall(0)
        dyn_conn=b_c.copy()
        conn=b_c.copy()
        G=nx.DiGraph()
        G.add_node(i)
        
        for j in range(len(self.BSpListBt)):
            # connectivity conditions
            if (not j==i) and self.GInBListBt[i][j]==0 and (self.BReacSpListBt[i] & self.BProdSpListBt[j]).any():
                dyn_conn[j]=1
                if not G.has_node(j):
                    G.add_node(j)
                
                if not G.has_edge(i,j):
                    G.add_edge(i,j,w=1)
                else:                        
                    G.edges[i, j]["w"]+=1
            elif (not j==i) and self.GInBListBt[i][j]==0 and (self.BProdSpListBt[i] & self.BReacSpListBt[j]).any():
                dyn_conn[j]=1
                if not G.has_node(j):
                    G.add_node(j)
                
                if not G.has_edge(j,i):
                    G.add_edge(j,i,w=1)
                else:                        
                    G.edges[j, i]["w"]+=1
            # Hasse condition (all connected)
            if (not j==i) and self.GInBListBt[i][j]==0:
                conn[j]=1
                
        return ([conn,dyn_conn,G])
    
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
    def setGeneratorsMp(self,threads=-1,verbose=False):
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
        self.RpInWhichGArray=-np.ones(self.MpDf.shape[1])
        self.RpInWhichGArray=self.RpInWhichGArray.tolist()
        # list of basic sets
        self.BSpListBt=[]
        
        # number of steps
        st=0
        # creating equivalence classes for each reaction that generate the same
        # closure and so each basic can be created
        
        for i in range(len(self.RpInWhichGArray)):
            if verbose:
                print(i+1," equivalance classes of ",len(self.RpInWhichGArray))
            
            st+=1
            # reaction already assigned to equivalnce class
            if self.RpInWhichGArray[i]>=0:
               continue
           
            self.RpInWhichGArray[i]=xeqc
            self.BSpListBt.append(c_reac[i])
            
            for j in range(i,self.MpDf.shape[1]):
                st+=1
                if c_reac[i]==c_reac[j]:
                    self.RpInWhichGArray[j]=xeqc
                    
            xeqc+=1
        
        # assingnig species contained each generator
        self.GSpListBt=[]
        
        b_sp=bt(self.MpDf.shape[0])
        b_sp.setall(0)
        
        for i in range(xeqc):
            self.GSpListBt.append(b_sp.copy())
            
        for i in range(self.MpDf.shape[1]):
            self.GSpListBt[self.RpInWhichGArray[i]]|= self.ReacListBt[i]|self.ProdListBt[i] 
        
        # assigning reaction contained each generator
        self.GRpListBt=[]
        b_r=bt(self.MpDf.shape[1])
        b_r.setall(0)
        
        for i in range(xeqc):
            for j in np.where(np.array(self.RpInWhichGArray)==i)[0]:
                b_r[j]=1   
            
            self.GRpListBt.append(b_r.copy())
            b_r.setall(0)
            
        
        if verbose:
            ver=100
        else: 
            ver=0
        
        parallel = Parallel( n_jobs = threads , require='sharedmem',verbose=ver)
        
        element_list=parallel( delayed( self.addGElements) ( i ) for i in range(len(self.GRpListBt)))
        
        self.GInBListBt=[]
        self.BRpListBt=[]
        self.BReacSpListBt=[]
        self.BProdSpListBt=[]
        self.BStoichioPositiveSpListBt=[]
        self.BStoichioNegativeSpListBt=[]
        
        for i in range(len(self.GRpListBt)):
            
            self.GInBListBt.append(element_list[i][0])
            self.BRpListBt.append(element_list[i][1])
            self.BReacSpListBt.append(element_list[i][2])
            self.BProdSpListBt.append(element_list[i][3])
            self.BStoichioPositiveSpListBt.append(element_list[i][4])
            self.BStoichioNegativeSpListBt.append(element_list[i][5])
        
        b_c=bt(len(self.BSpListBt))
        b_c.setall(0)

        conn_list=parallel( delayed( self.connGGen) ( i ) for i in range(len(self.GRpListBt)))
        
        self.NotContainedBListBt=[]
        self.ConnectedBListBt=[]
        self.GConnNx=nx.DiGraph()
        
        for i in range(len(self.GRpListBt)):
            self.NotContainedBListBt.append(conn_list[i][0])
            self.ConnectedBListBt.append(conn_list[i][1])
            Gcomp = nx.compose(self.GConnNx, conn_list[i][2])
            edge_data = {e: self.GConnNx.edges[e]['w'] + conn_list[i][2].edges[e]['w'] for e in self.GConnNx.edges & conn_list[i][2].edges}
            nx.set_edge_attributes(Gcomp, edge_data, 'w')
            self.GConnNx=Gcomp.copy()
                
        return
    
    def MgenByInd(self,i):
        '''
        

        Parameters
        ----------
        i : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
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
        
        return(v)
        
    # Minimal generators generation function, to be started once the basic sets
    # have been generated by the setGenerators() function.  It returns a list
    # (MgenListListSpArray) of bitarray list of species that correspond to the sets so that 
    # by means of the closure they generate the basic molecule.
    def setMgenMp(self,threads=-1,verbose=False):
        if verbose:
            ver=100
        else:
            ver=0
        if not hasattr(self, 'GInBListBt'):
            print("The basic sets have not been initialized, please run the setGenerators() function.")
            return 
        self.MGenStepInt=0
        parallel = Parallel( n_jobs = threads , require='sharedmem',verbose=ver)
        self.MgenListListSpArray=parallel( delayed( self.MgenByInd) ( i ) for i in range(len(self.GRpListBt)))
    
    # Synergy generation function. it takes as input a minimum generator (sp)
    # and (pi) the index of the objective . 
    # The function verifies which combinations of generators can 
    # generate such a generator.        
    def genSynMp(self,sp,pi,ReacSyn,ProdSyn):

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
            self.recursiveGenSynMp(p,i,sp,xp,pi,ReacSyn,ProdSyn)
        
        return
    
    # Recursive synergy generation function, requires as inputs (p) the 
    # existing generators to combine, o the next level to add to the scan, 
    # xp list of indexes of the corresponding generators, (sp) minimum 
    # generator to reach and pi the index of the objective generator. 
    # Function recursively explores the possible 
    # combinations to reach a synergy. If this is reached, the branch 
    # to be explored will be cut. The synergies are stored in a list (SynReacListGBt).
    def recursiveGenSynMp(self,p,o,sp,xp,pi,ReacSyn,ProdSyn):
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
                    u|=self.GSpListBt[xp[j]]
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
                c_syn[xp[i]]=1
            
            # if not c_syn in ReacSyn:
            ReacSyn.append(c_syn)
            op=bt(len(self.GInBListBt))
            op.setall(0)
            op[pi]=1
            ProdSyn.append(op)

                
            p[o]=0
            return
        
        # if not recursion continue
        for i in range(o+1,len(xp)):
            self.recursiveGenSynMp(p,i,sp,xp,pi,ReacSyn,ProdSyn)
     
        p[o]=0
    
    def SynByInd(self,i):
        '''
        

        Parameters
        ----------
        i : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        ReacSyn=[]
        ProdSyn=[]
        for j in self.MgenListListSpArray[i]:
            if j.count()>0:
                self.genSynMp(j,i,ReacSyn,ProdSyn)
                # ReacSyn+=list(filter(lambda x: not x in ReacSyn, TempReacSyn))
                # ProdSyn+=list(filter(lambda x: not x in ProdSyn, TempProdSyn))
                
            else:
                for k in range(len(self.MgenListListSpArray)):
                    if k!=i:
                        op=bt(len(self.GInBListBt))
                        op.setall(0)
                        op[k]=1
                        ReacSyn.append(op.copy())
                        op.setall(0)
                        op[i]=1
                        ProdSyn.append(op.copy())
        return[ReacSyn,ProdSyn]

    # Function that generates all the synergies from the minimum 
    # generators. This is achieved through the use of the function genSyn()
    # and the recursive function recursiveGenSyn(). The output consists of list SynReacListGBt 
    # which contains all the synergies and list (SynProdListGBt) which contains 
    # all the triggered generators.   
    def setSynMp(self,threads=-1,verbose=False):
        '''
        

        Parameters
        ----------
        threads : TYPE, optional
            DESCRIPTION. The default is -1.

        Returns
        -------
        None.

        '''
        if not hasattr(self, 'MgenListListSpArray'):
            print("The minimal genetators have not been initialized, please run the setMgen(() function.")
            return         
        # List of posible synergies
        if verbose:
            ver=100
        else:
            ver=0
            
        parallel = Parallel( n_jobs = threads ,verbose=ver)
        all_syn=parallel( delayed( self.SynByInd) ( i ) for i in range(len(self.MgenListListSpArray)))
        
        self.SynReacListGBt=[]
        self.SynProdListGBt=[]
        
        for i in all_syn:
            self.SynReacListGBt+=i[0]
            self.SynProdListGBt+=i[1]
    

    def genSynClos(self,i):
        '''
        

        Parameters
        ----------
        Gset : bitarray
            Proto-synergy reactive part bitarray.

        Returns
        -------
        c : bitarray
            Closure of proto-synergy as bitarray of generators.

        '''
        c=self.getGBtInSpBt(self.getClosureFromSp(self.getSpBtInGBt(self.SynProdListGBt[i]|self.SynReacListGBt[i]),bt_type=True))
        return c
        
    def setSynEcMp(self,threads=-1,verbose=False):
        '''
        

        Parameters
        ----------
        threads : TYPE, optional
            DESCRIPTION. The default is -1.

        Returns
        -------
        None.

        '''
        if verbose:
            ver=100
        else:
            ver=0
            
        # creating all closed set of closure of reactant synegric part 
        parallel = Parallel( n_jobs = threads , require='sharedmem',verbose=ver)
        syn_c_reac=parallel( delayed( self.genSynClos) ( i ) for i in range(len(self.SynReacListGBt)))
        
        # class equivalence Minimal synergistic clousrues dictonary
        self.MinSynCLDict={}
        
        # Creating the dictonary for the proto-synergies reactive parts
        for i in range(len(self.SynReacListGBt)):
            if verbose:
                print(i+1," proto-synergies of ",len(self.SynReacListGBt))
            
            # reaction already assigned to equivalnce class
            if fbt(syn_c_reac[i]) in self.MinSynCLDict.keys():
               self.MinSynCLDict[fbt(syn_c_reac[i])].append(self.SynReacListGBt[i])
            else:
                self.MinSynCLDict[fbt(syn_c_reac[i])]=[self.SynReacListGBt[i]]
            
    def initSynStr(self,i,G):
        '''
        

        Parameters
        ----------
        i : TYPE
            DESCRIPTION.
        G : TYPE
            DESCRIPTION.

        Returns
        -------
        
        '''
        is_ssm=self.isSsmFromSp(self.BSpListBt[i])
        if is_ssm:
            is_org=self.isSmFromSp(self.BSpListBt[i])
            G.add_node(fbt(self.GInBListBt[i]),level=self.GInBListBt[i].count(),
                       sp=self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])],
                       is_ssm=is_ssm,is_org=is_org,is_basic=True,basic_id=i)
        else:    
            G.add_node(fbt(self.GInBListBt[i]),level=self.GInBListBt[i].count(),
                       sp=self.SpIdStrArray[self.getIndArrayFromBt(self.BSpListBt[i])],
                       is_ssm=is_ssm,is_org=False,is_basic=True,basic_id=i)
     
        
    def getClosFromG(self,GSet):
        '''
        

        Parameters
        ----------
        GSet : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        
        MinSynCL = list(map(lambda x: x[0], filter(lambda x: any(list(map(lambda y: y==(GSet & y), x[1]))) and (x[0] != GSet & x[0]), self.MinSynCLDict.items())))
        
        if len(MinSynCL) > 0:
           MinSynCL = list(filter(lambda x: x.count() == max(list(map(count),MinSynCL))))         
           return MinSynCL[0]
        else:
            return GSet
        
    def synStrPerNode(self,j):
        G=nx.MultiGraph()
        # Generating closures whit connected basics sets for each set in level i
        for k in self.getIndArrayFromBt(self.getGBtNotInBBt(bt(j))): 
            # Closure result
            add_gen=bt(j)
            add_gen.setall(0)
            add_gen[k]=1
            cr_a=fbt(self.getClosFromG(add_gen|j))
            cr_sp=self.getSpBtInGBt(cr_a)
            
            # node is added if is not in structrue
            if not (cr_a in G):
                is_ssm=self.isSsmFromSp(cr_sp)
                if is_ssm:
                    is_org=self.isSmFromSp(cr_sp)
                    G.add_node(cr_a,level=cr_a.count(),
                               sp=self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)],
                               is_ssm=is_ssm,is_org=is_org,is_basic=False)
                else:
                    G.add_node(cr_a,level=cr_a.count(),
                               sp=self.SpIdStrArray[self.getIndArrayFromBt(cr_sp)],
                               is_ssm=is_ssm,is_org=False,is_basic=False)

            # Adding edges corresponding to the colsure, and verifing if is a sinergy:
            if cr_a.count() > (bt(j)|self.GInBListBt[k]).count():
               G.add_edge(j,cr_a,key=fbt(self.GInBListBt[k]),syn=True,added_basic=k)
            else:
               G.add_edge(j,cr_a,key=fbt(self.GInBListBt[k]),syn=False,added_basic=k)
           # if cr_a.count() > (bt(j)|self.GInBListBt[k]).count():
               #    G.add_edge(j,cr_a,key=fbt(self.GInBListBt[k]),syn=True,added_basic=k)
               # else:
               #    G.add_edge(j,cr_a,key=fbt(self.GInBListBt[k]),syn=False,added_basic=k)
        return G
     
        
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
    def setSynStrMp(self,partial_save=None,threads=-1,verbose=False):
        if verbose:
            ver=100
        else:
            ver=0
        if not hasattr(self, 'GInBListBt'):
            print("The basic sets have not been initialized, please run the setGenerators() function.")
            return 
        # Initialization of the synergistic structure as a multigraph object
        G = nx.MultiDiGraph()
        parallel = Parallel( n_jobs = threads , require='sharedmem',verbose=ver)
        parallel( delayed( self.initSynStr) ( i, G ) for i in range(len(self.BSpListBt)))
                
        # Generation of the multigraph of the synergistic structure by set level (contained basics)
        for i in range(len(self.BSpListBt)):
            
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
            
            # Generating closures whit connected basics sets for each set in level i
            parallel = Parallel( n_jobs = threads , require='sharedmem',verbose=ver)
            G_list=parallel( delayed( self.synStrPerNode) ( i ) for i in nodes)
            for i in G_list:
                G = nx.compose(G, i)
            
            if not partial_save is None:
                self.saveToPkl(partial_save)
        
        self.SynStrNx=G
