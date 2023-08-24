#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 14:07:20 2022

@author: pmaldona

Reaction Network Decomposition Class 
"""
from .RNIRG import RNIRG
import numpy as np
from bitarray import bitarray as bt
from bitarray import frozenbitarray as fbt
from scipy.optimize import linprog
import networkx as nx
from itertools import chain, combinations
from typing import List, Any, Iterable
from functools import reduce
from bitarray.util import subset
from pulp import *
import re
import random
from pyvis.network import Network
from colorsys import hsv_to_rgb
import pypoman as ph

class RNDS(RNIRG):

            
    def getOpOrgProcess(self,sp_set,opsp_set=None,pr=None):
        '''
        

        Parameters
        ----------
        sp_set : Numpy arrary or Bitarray
            Species Id array or bitarray of species of an organization.
        opsp_set : Numpy arrary or Bitarray, 
            Species Id Array or species bitarray of species that must be 
            overproduced. The default is None.
        pr : Numpy array ot bitarray, optional
            indexes of active reactions. The default is None.

        Returns
        -------
        sp : list
            A feasible procees vector for a organization that also overproduce at least opsp_set.

        '''
        # Checks if input is bitarray, if it is not, it make the 
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
        
        if opsp_set is None:
            opts=sp.copy()
        else:
            if not (isinstance(opsp_set,bt)):
                opts=bt(self.MpDf.shape[0])
                opts.setall(0)
                
                for i in opsp_set:
                    if i in self.MpDf.index.values:
                        ind=self.MpDf.index.get_loc(i)
                        opts[ind]=1
            else:
                opts=opsp_set
        
        if pr is None:
            v=bt(self.MpDf.shape[1])
            v.setall(1)    
            
        elif not (isinstance(pr,bt)):
            v=bt(self.MpDf.shape[1])
            v.setall(0)
            
            for i in pr:
                v[i]=1
        else:
            v=pr
            
        # If it's only one species and self mantained, the reuturns the specie
        # itself
        if sp.count()==1 and self.isSmFromSp(sp):
            return sp
        
        Ns=sp.count() # number of species (rows)
        nsp=list(set(range(len(sp)))-set(self.getIndArrayFromBt(sp))) #species no present
        
        rc =[] # creating variable of available reaction 
        for j in self.getIndArrayFromBt(v):
            if (all(self.MpDf.iloc[nsp,j]==0) and all(self.MrDf.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
        
        if len(rc)==0:
            sp.setall(0)
            return sp
        
        
        # stoichiometric matrix of rn with prepended destruction reactions
        S=self.MpDf.iloc[self.getIndArrayFromBt(sp),rc]-self.MrDf.iloc[self.getIndArrayFromBt(sp),rc]
        sp_names=S.index
        S=-S.to_numpy()
        S=S.tolist()
        
        # opverprodcible index respect a reduced stochimetric matrix
        op_ind=np.empty(opts.count())
        k=0
        for i in opts.search(1):
            op_ind[k]=np.where(self.SpIdStrArray[i] == sp_names)[0]
            k+=1
        op_ind=op_ind.astype(int)
        
        
        # flow vector constraint: production of every species = 0
        f=np.zeros(Ns) #norm of porcess vector for minimization 
        f[op_ind]=-1 # condition for the overproducied species
        f=f.tolist()
        
        # lower bound for all reaction, so all have to take place
        bounds=[]
        for i in range(len(rc)):
            bounds.append((1,None))
           
        # cost 0 for every reaction
        cost=np.ones(len(rc)) #norm of porcess vector for minimization 
        cost=cost.tolist()
        
        # Resloving the Lp under the organization condition 
        res = linprog(cost, A_ub=S, b_ub=f, bounds=bounds, method='highs')
        
        # creating result vector
        v=np.zeros(self.MpDf.shape[1])
        if res.status==0:
            k=0
            for i in rc:
               v[i]=res.x[k]
               k+=1
            return(v)
        
        else:
            return(None)

    
    def veriOpSpBt(self,sp_set,opsp_set,pr=None,force_org=False):     
        '''
        

        Parameters
        ----------
        sp_set : bitarray
            Set of species to check.
        opsp_set : bitarray
            Set of overproduced whitin sp_set to check..
        pr : list of integers or numpy array of integers, optional
            Indexes of active reactions. The default is None.
        force_org: bool
            Forces the activacion all reactions (rates greater that zero) of pr

        Returns
        -------
        bool
            True if opsp are overproduced in the set sp_set, and reaction pr. Else False.
        numpy array
            If latter condition is True, returns related process vector. If False, return vector of zeros.

        '''
        # Checks if input is bitarray, if it is not, it make the 
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
        
        if opsp_set is None:
            opts=sp.copy()
        else:
            if not (isinstance(opsp_set,bt)):
                opts=bt(self.MpDf.shape[0])
                opts.setall(0)
                
                for i in opsp_set:
                    if i in self.MpDf.index.values:
                        ind=self.MpDf.index.get_loc(i)
                        opts[ind]=1
            else:
                opts=opsp_set
        
        if pr is None:
            v=bt(self.MpDf.shape[1])
            v.setall(1)    
            
        elif not (isinstance(pr,bt)):
            v=bt(self.MpDf.shape[1])
            v.setall(0)
            
            for i in pr:
                v[i]=1
        else:
            v=pr
        
        opsp=self.SpIdStrArray[opsp_set.search(1)]
        
        nsp=list(set(range(len(sp)))-set(self.getIndArrayFromBt(sp))) #species no present
           
        rc =[] # creating variable of available reaction 
        rsp=set()
        for j in self.getIndArrayFromBt(v):
            if (all(self.MpDf.iloc[nsp,j]==0) and all(self.MrDf.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
                rsp=rsp.union(np.where(self.MpDf.iloc[:,j]!=0)[0]).union(np.where(self.MrDf.iloc[:,j]!=0)[0])
        
        rsp=list(rsp)
        rsp=self.SpIdStrArray[rsp]
# 
        if len(rc)==0:
            return False, np.zeros(self.MpDf.shape[1])
        
     
        else:   
               
            # stoichiometric matrix of rn with prepended destruction reactions
            S=self.MpDf.iloc[self.getIndArrayFromBt(sp),rc]-self.MrDf.iloc[self.getIndArrayFromBt(sp),rc]
            # S=np.column_stack((-np.identity(Ns),S))
            # S=S.tolist()
               
            # flow vector constraint: production of every species = 0
            # f=np.zeros(Ns) #norm of porcess vector for minimization 
            # f=f.tolist()
               
            # Name and type of problem
            prob = LpProblem("LP_Org_test", LpMinimize)
               
            # Definition of problem warabiles
            if force_org:
                x = LpVariable.dicts("process", S.iloc[0].index, lowBound=1, cat='Continuous')
            else:
                x = LpVariable.dicts("process", S.iloc[0].index, lowBound=0, cat='Continuous')
          
            
            # We define the objective function
            prob += lpSum([x[i] for i in x])
            for j in S.index:
                    if j in opsp:
                        prob +=lpSum([S.loc[j][i] * x[i] for i in S.loc[j].index]) >= 1
                    elif j in rsp:
                        prob +=lpSum([S.loc[j][i] * x[i] for i in S.loc[j].index]) == 0            

                  
            # Resolver el problema
            prob.solve(solver=GLPK(msg=False))
               
            v_result=np.zeros(self.MpDf.shape[1])
            for v in prob.variables():
                if not ("dummy" in v.name):
                    ind=int(re.search(r'\d+', v.name).group())
                    v_result[ind]=v.varValue
                
            if LpStatus[prob.status]=='Optimal':
                return True, v_result
            else:
                return False, np.zeros(self.MpDf.shape[1])
            
            
    # Verifies which species of a reaction network (not necessarily an organization) 
    # are overproducible. Inputs are species sp_set and process vector pr, 
    # returns a list of overproducible species
    def getallOpSpBt(self,sp_set,pr=None,):     
        # Checks if input is bitarray, if it is not, it make the 
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
        
        
        if pr is None:
            v=bt(self.MpDf.shape[1])
            v.setall(1)    
            
        elif not (isinstance(pr,bt)):
            v=bt(self.MpDf.shape[1])
            v.setall(0)
            
            for i in pr:
                v[i]=1
        else:
            v=pr
            
        # If it's only one species and self mantained, the reuturns the specie
        # itself
        if sp.count()==1 and self.isSmFromSp(sp):
            return sp
        
        Ns=sp.count() # number of species (rows)
        nsp=list(set(range(len(sp)))-set(self.getIndArrayFromBt(sp))) #species no present
        
        rc =[] # creating variable of available reaction 
        for j in self.getIndArrayFromBt(v):
            if (all(self.MpDf.iloc[nsp,j]==0) and all(self.MrDf.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
        
        if len(rc)==0:
            sp.setall(0)
            return sp
        
        
        # stoichiometric matrix of rn with prepended destruction reactions
        S=self.MpDf.iloc[self.getIndArrayFromBt(sp),rc]-self.MrDf.iloc[self.getIndArrayFromBt(sp),rc]
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
                opsp[self.getIndArrayFromBt(sp)[i]]=1
                
        return(opsp)
    
    # Verifies which species of a reaction network (not necessarily an organization) 
    # are overproducible. Inputs are species sp_set and process vector pr, 
    # returns a list of overproducible species
    def getSpNeededToOrg(self,sp_set,pr=None,destruct=False):     
        # Checks if input is bitarray, if it is not, it make the 
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
        
        if pr is None:
            v=bt(self.MpDf.shape[1])
            v.setall(1)    
            
        elif not (isinstance(pr,bt)):
            v=bt(self.MpDf.shape[1])
            v.setall(0)
            
            for i in pr:
                v[i]=1
        else:
            v=pr
            
        # If it's only one species and self mantained, the reuturns the specie
        # itself
        if sp.count()==1 and self.isSmFromSp(sp):
            return sp
        
        Ns=sp.count() # number of species (rows)
        nsp=list(set(range(len(sp)))-set(self.getIndArrayFromBt(sp))) #species no present
        
        rc =[] # creating variable of available reaction 
        for j in self.getIndArrayFromBt(v):
            if (all(self.MpDf.iloc[nsp,j]==0) and all(self.MrDf.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
       
        # stoichiometric matrix of rn with prepended creation and destruction reactions 
        S=self.MpDf.iloc[self.getIndArrayFromBt(sp),rc]-self.MrDf.iloc[self.getIndArrayFromBt(sp),rc]
        S=S.to_numpy()
        S=np.column_stack((np.identity(Ns),np.column_stack((-np.identity(Ns),S))))
        #print("S:",S)
        S=S.tolist()
       
        
        # production of every species = 0, always possible because of additional creation and destruction reactions
        f=np.zeros(Ns) 
        f=f.tolist()
        #print("f:",f)
        # original reactions with rate>=1 (any positive), prepended reactions with rate>=0
        h=-np.ones(2*Ns+len(rc))
        h[0:2*Ns]=0
        h=h.tolist()
        #print("h:",h)
        
        # cost 0 for every reaction
        cost=np.zeros(2*Ns+len(rc)) #norm of porcess vector for minimization
        cost[0:Ns+Ns*destruct]=1
        cost=cost.tolist()
        #print("cost:",cost)
        
        G=-np.identity(2*Ns+len(rc))
        #print("G:",G)
        G=G.tolist()
        
        res = linprog(cost, A_eq=S, b_eq=f, A_ub=G, b_ub=h, method='highs')
        # The unknown to be solved is the process column vector v = [vc,vd,vo] (creation/destruction/original network).
        # The result for v is stored in the variable rn.linp.r$X (vc is X[1:Ns], vd is X[Ns+(1:Ns)]).
        # The equations are S v = f (with f=0), G v >= h (i.e. [vc,vd,vo] >= [0,0,1] because G is the identity matrix).
        # Only the original network reactions are constrained to be strictly positive (here arbitrarily vo >= 1).
        # The cost is  Cost . vc  (dot product) or  Cost . vd  when destruction is penalized instead of creation.
        # The linear programming ideal result is every prepended creation reaction with rate 0 (minimal Cost 0).
        # It implies that there is at least one original process vector that can sustain every species (and even increase
        # them if destruction reactions have positive rate). In that case the reaction network is a proper organization.
        # If not it can be converted into one by adding the extra inflow of species with creation rate > 0.
        # Returns the indexes of species needed as extra inflow and overproduced species in the particular found solution.        
        
        # return(dict(inflow=np.where(res.x[0:Ns]>0)[0],opsp=np.where(res.x[Ns+1:2*Ns]>0)[0],v=res.x[0:Ns]))
        return (np.where(res.x[0:Ns]>0)[0],np.array(rc))


    # Generates a base of overproduced species of a reaction network 
    # (not necessarily an organization). Inputs are species sp_set and 
    # process vector pr, returns a list of minium overproducible species sets.    
    def getOpBaseBtList(self,sp_set,pr=None):     
        # Checks if input is bitarray, if it is not, it make the 
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
        
        if pr is None:
            v=bt(self.MpDf.shape[1])
            v.setall(1)    
            
        elif not (isinstance(pr,bt)):
            v=bt(self.MpDf.shape[1])
            v.setall(0)
            
            for i in pr:
                v[i]=1
        else:
            v=pr
            
            
        # If it's only one species and self mantained, the reuturns the specie
        # itself
        if sp.count()==1 and self.isSmFromSp(sp):
            return {"op_b":[sp]}
        
        Ns=sp.count() # number of species (rows)
        nsp=list(set(range(len(sp)))-set(self.getIndArrayFromBt(sp))) #species no present
        
        rc =[] # creating variable of available reaction 
        for j in self.getIndArrayFromBt(v):
            if (all(self.MpDf.iloc[nsp,j]==0) and all(self.MrDf.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
       
        # stoichiometric matrix of rn with prepended destruction reactions
        S=self.MpDf.iloc[self.getIndArrayFromBt(sp),rc]-self.MrDf.iloc[self.getIndArrayFromBt(sp),rc]
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
                        opsp[self.getIndArrayFromBt(sp)[i]]=1
                
                if not (opsp in op_b):
                    op_b.append(opsp)
                    
                    # adding corresponding proces to the base
                    # vi=np.zeros(self.MpDf.shape[1])
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
    def getDcomArray(self,sp_set,opsp_set,pr=None):
        
        # Checks if input is bitarray, if it is not, it make the 
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
            
        if not (isinstance(opsp_set,bt)):
            opsp=bt(self.MpDf.shape[0])
            opsp.setall(0)
            
            for i in opsp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    opsp[ind]=1
        else:
            opsp=opsp_set
      
        if pr is None:
            v=bt(self.MpDf.shape[1])
            v.setall(1)    
            
        elif not (isinstance(pr,bt)):
            v=bt(self.MpDf.shape[1])
            v.setall(0)
            
            for i in np.where(pr>0)[0]:
                v[i]=1
        else:
            v=pr
            
        # generetinf the decompotion vector
        dcom=np.zeros(len(sp))
        
        # # generating overproducible species
        # opsp = self.over_prod(sp,v)
        dcom[self.getIndArrayFromBt(opsp)]=-1
        
        # If it's only one species and self mantained, the reuturns the specie
        # itself
        if sp.count()==1 and self.isSmFromSp(sp):
            return dcom
        
        
        # generating matrices to find catalytic and non catalytic species
        sp_ind=self.getIndArrayFromBt(sp)
        c_m = self.MrDf.iloc[sp_ind,:].copy()
        nc_m = c_m.copy()
        
        #not present specues as index
        nsp=sp.search(0)
        
        rsp=set()
        rc=[]
        for i in self.getIndArrayFromBt(v):
            # Matrix for identifaction of catalysts species identifications
            c_m.iloc[:,i] = ((self.MpDf.iloc[sp_ind,i]!=0) & (self.MrDf.iloc[sp_ind,i]!=0)) & (self.MpDf.iloc[sp_ind,i]==self.MrDf.iloc[sp_ind,i])
            # Matrix for identifaction of non-catalysts species identifications
            nc_m.iloc[:,i] = (self.MpDf.iloc[sp_ind,i])!=(self.MrDf.iloc[sp_ind,i])
            
            if (all(self.MpDf.iloc[nsp,i]==0) and all(self.MrDf.iloc[nsp,i]==0)):
                rc.append(i) # selecting reactions that can be trigger with available species
                rsp=rsp.union(np.where(self.MpDf.iloc[:,i]!=0)[0]).union(np.where(self.MrDf.iloc[:,i]!=0)[0])
        
        # Removing reaction with zero rate
        c_m=c_m.iloc[:,self.getIndArrayFromBt(v)]
        nc_m=nc_m.iloc[:,self.getIndArrayFromBt(v)]     
         
        # finding catalytic species
        csp=sp.copy()
        csp.setall(0)
        j=0 
        for i in ((c_m.sum(axis=1)!=0) & (nc_m.sum(axis=1)==0)):
            if i:
                csp[sp_ind[j]]=1
            j+=1
        
        # non-reactive species
        nrsp_ind=list(set(self.getIndArrayFromBt(sp))-rsp)

        nrsp=sp.copy()
        nrsp.setall(0)
        
        # as index vector
        for i in nrsp_ind:
            nrsp[i]=1
            # setting as -3 non reactive species
            dcom[i]=-3
        
        # findinng all other species that aren't catalysers o opverproducble
        fsp=sp.copy() & ~(csp.copy() | opsp.copy() | nrsp.copy())
        
        # if there aren't reactive 
        if(len(rc)==0):
            return
        
        # creating fragile cycles
        if fsp.count()>0:            
            # subselecting matrix for fragile cycles species
            m=self.MpDf.iloc[self.getIndArrayFromBt(fsp),rc]-self.MrDf.iloc[self.getIndArrayFromBt(fsp),rc]
            
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
                dcom[self.getIndArrayFromBt(fsp)[i]]=fsp_v[i]
                
        # identifying catalyst species
        dcom[self.getIndArrayFromBt(csp)]=-2
        return dcom
                    
                                       
    # Function that verifiys if a decomposition input is semi-self-
    # mantained or not.
    def getSsmDcomArray(self,dcom):
        
        #classification of overproduced species
        opsp=np.where(dcom==-1)[0]

        nsp=np.where(dcom==0)[0] #no present species
        psp=np.where(dcom!=0)[0] #present species
        
        rc =[] # creating variable of available reaction 
        for j in range(self.MpDf.shape[1]):
            if (all(self.MpDf.iloc[nsp,j]==0) and all(self.MrDf.iloc[nsp,j]==0)):
                rc.append(j) # selecting reactions that can be trigger with available species
        # if there aren't reactive 
        if(len(rc)==0):
            return
        
        # Stochimetric matrix whit available species and reactions
        S=-self.MpDf.iloc[psp,rc]+self.MrDf.iloc[psp,rc]
        S=S.values.tolist()
        
        # condition of overproduction for overproduced species
        f=np.zeros(len(psp))
        f[np.isin(psp,opsp)]=-1
        f=f.tolist()
        
        nfcsp=list(set(nsp) | set(opsp)) #non present and opverproduced species, 
        fcrc =[] # creating variable of available reaction for fragile cycles 
        for j in range(self.MpDf.shape[0]):
            if (all(self.MpDf.iloc[nfcsp,j]==0) and all(self.MrDf.iloc[nfcsp,j]==0)):
                fcrc.append(j) # selecting reactions that can be trigger with available species
        # if there aren't reactive 
        # Matrix generation for the equality so fragile cycles species have production 0
        S_eq=(self.MpDf-self.MrDf).to_numpy()
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
    
    
    def getallPermutation(self,sp_set,level):
        '''
        

        Parameters
        ----------
        sp_set : TYPE
            DESCRIPTION.
        level : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        all_comb=combinations(self.getIndArrayFromBt(sp_set), level)
        return list(map(lambda x: fbt(self.getBtFromIndArray(x,len(sp_set))),all_comb))
    
    def addfOpfromList(self,G,sp,oplist,pr,force_org=False):
        '''
        

        Parameters
        ----------
        G : TYPE
            DESCRIPTION.
        sp : TYPE
            DESCRIPTION.
        oplist : TYPE
            DESCRIPTION.
        pr : TYPE
            DESCRIPTION.
        force_org : bool
            Forces the activacion all reactions (rates greater that zero) of pr

        Returns
        -------
        None.

        '''
        for i in oplist:
            # Verifying if propused set of overproducibles are overproduced
            verifi , process = self.veriOpSpBt(sp_set=sp,opsp_set=bt(i),pr=pr,force_org=force_org)
            
            # Estimating whether lower level sets exist
            try:
                low_level=list(map(lambda x: bt(x[0]),filter(lambda x: x[1]['level']<=i.count()-1, G.nodes(data=True))))
                inf_level_sp=reduce(lambda x,y : x|y,low_level)
            except:
                inf_level_sp=bt(i)
                inf_level_sp.setall(0)
            
            # Adding new nodes
            if verifi and (not subset(i, inf_level_sp)):
                dcom=self.getDcomArray(sp, i, pr)
                G.add_node(i, level=i.count(), process=process,decomposition=dcom)
              
                # Adding new edges to hasse diagram
                for m in low_level:    
                    if m & i == m:  
                        G.add_edge(fbt(m),fbt(i))
            
              
    def genOpBase(self,sp_set=None,pr=None,force_org=False):
        '''
        

        Parameters
        ----------
        sp_set : TYPE, optional
            DESCRIPTION. The default is None.
        pr : TYPE, optional
            DESCRIPTION. The default is None.
        force_org : bool
            Forces the activacion all reactions (rates greater that zero) of pr
        Returns
        -------
        G : TYPE
            DESCRIPTION.

        '''
        
        if sp_set is None:
            sp=bt(self.MpDf.shape[0])
            sp.setall(1)
            
        elif not (isinstance(sp_set,bt)):
            sp=bt(self.MpDf.shape[0])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.MpDf.index.values:
                    ind=self.MpDf.index.get_loc(i)
                    sp[ind]=1
        else:
            sp=sp_set
        
        if pr is None:
            v=bt(self.MpDf.shape[1])
            v.setall(1)    
            
        elif not (isinstance(pr,bt)):
            v=bt(self.MpDf.shape[1])
            v.setall(0)
            
            for i in pr:
                v[i]=1
        else:
            v=pr
        
        opsp=self.getallOpSpBt(sp,v)
        G = nx.DiGraph() # lattice overproducible structre as network graph
        # We want to generate al possible overproducible species combiantion 
        # that are ovreproduced but not comination of existing overprducible species
        for level in range(1,opsp.count()+1):
            
            # try:
            #     sup_level_sp=reduce(lambda x,y : x|y,map(lambda x: bt(x[0]),filter(lambda x: x[1]['level']>level, G.nodes(data=True))))
            # except:
            #     break
            try:
                # we get the union of all inferior species if it is possible (case we are not at level=1)
                inf_level_sp=reduce(lambda x,y : x|y,map(lambda x: bt(x[0]),filter(lambda x: x[1]['level']<=level-1, G.nodes(data=True))))
                pos_comb=self.getallPermutation(inf_level_sp, level)
            except:
                pos_comb=[]
                
            all_comb=self.getallPermutation(opsp, level) #we obtain all posible combinations
            new_comb=list(filter(lambda x: x not in pos_comb, all_comb))
            # new_comb=list(filter(lambda x: x not in list(G.nodes), new_comb))
            if len(new_comb)>0: 
                self.addfOpfromList(G, sp, new_comb, pr,force_org=force_org)                        
            
            else:
                # No more cases to explore
                break
            
        return G
         
                    
            
    # Function that generates a Hasse diagram from a set of species with all 
    # the combinations of the overproduced base, obtaining all the 
    # decompositions. The input corresponds to a set of sp_set species which 
    # must be an organization and the output corresponds to a graph whose 
    # vertexs correspond to the different combinations of overproduced species. 
    # The edges corresponds to the join and meet of the overproduced species. 
    # Each vertex has an attribute that is the decomposition vector. 
    def getOpHasseNx(self,sp_set,force_org=False):
        # Checks if input is a bitarray, if it is not, it make the 
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
      
        # initialization of multigraph of opverproduced by inducing the base 
        op_hasse = self.genOpBase(sp_set, self.getTriggerableRpBtFromSp(sp_set),force_org=force_org)
        all_op=sp.copy()
        all_op.setall(0)
        
        op_b=[]
        
        # The nodes corresponding to the overproduced base.
        for i in op_hasse.nodes():
            op_b.append(bt(i))
            all_op|=bt(i)
        
        # print("total levels: ",all_op.count())  
        # Generation of the multigraph of overproduced hasse, by number of 
        # overproduced species
        for i in range(1,all_op.count()+1):
            level_nodes=list(map(lambda x: bt(x[0]),filter(lambda x: x[1]['level']==i, op_hasse.nodes(data=True))))
            for j in level_nodes:
                for k in op_b:
                    op_new= fbt(j|k)
                    # print("exploring: ",op_new)
                    if not (op_new in op_hasse.nodes()):
      
                        process=op_hasse.nodes(data="process")[fbt(j)]+op_hasse.nodes(data="process")[fbt(k)]
                        dcom=self.getDcomArray(sp,bt(op_new),pr=np.where(process>0)[0])
                        op_hasse.add_node(fbt(op_new),level=op_new.count(),
                                          decomposition=dcom,
                                          process=process)#,
         
                        if op_new.count()>0:
                            low_level_nodes = list(map(lambda x: x[0],filter(lambda x: x[1]['level']==op_new.count()-1, op_hasse.nodes(data=True))))
                            # [x for x,y in op_hasse.nodes(data=True) if y['level']==i.count()-1]
                            for m in low_level_nodes:    
                                if bt(m) & bt(op_new) == bt(m):  
                                    op_hasse.add_edge(m,fbt(op_new))
                 
        return op_hasse
    
    def getDecomDisplayPv(self,decom_array,process,x_size='500px',y_size='500px',notebook=False,cdn_resources='local',disp_non_act=True,sp_name=False):

        G = nx.MultiDiGraph()
        r_i=np.where(process>0)[0]
        
        exluded_colors=["#E0E0E0","#0080FF","#FF6666","#00FF80","#4B0082"]
        
        if sp_name:
            sp_id=self.SpNameStrArray.copy()
        else: 
            sp_id=self.SpIdStrArray.copy()
            
        for i in np.unique(decom_array):
            if i ==-3 and not disp_non_act:
                for j in np.where(decom_array==i)[0]:
                    G.add_node(self.SpIdStrArray[j], color = "#FF6666", label=sp_id[j], size=14, shape="dot")
            elif i == -2:
                for j in np.where(decom_array==i)[0]:
                    G.add_node(self.SpIdStrArray[j], color = "#0080FF", label=sp_id[j], size=14, shape="dot")
            elif i == -1:
                for j in np.where(decom_array==i)[0]:
                    G.add_node(self.SpIdStrArray[j], color = "#00FF80", label=sp_id[j], size=14, shape="dot")
            elif i == 0 and disp_non_act:
                for j in np.where(decom_array==i)[0]:
                    G.add_node(self.SpIdStrArray[j], color = "#E0E0E0", label=sp_id[j], size=14, shape="dot")
            elif i>0:                
                # generate a new color that is visually distinguishable from the existing colors
                while True:
                    hue = random.random()
                    saturation = random.uniform(0.5, 1.0)
                    brightness = random.uniform(0.5, 1.0)
                    r, g, b = hsv_to_rgb(hue, saturation, brightness)
                    color = "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))
                    if color not in exluded_colors:
                        break
                    
                for j in np.where(decom_array==i)[0]:
                    G.add_node(self.SpIdStrArray[j], color = color, label=sp_id[j], size=14, shape="dot")
        
        for i in range(self.MpDf.shape[1]):
            if i in r_i:
                if process[i]==1:
                    G.add_node("r"+str(i), color = "#4B0082", label="r"+str(i), shape="square", size=7,title="")
                else:
                    G.add_node("r"+str(i), color = "#4B0082", label="r"+str(i), shape="square", size=7,title=str(process[i]))
                    
                for j in self.getIndArrayFromBt(self.ReacListBt[i]):
                    st_value=self.MrDf.iloc[j,i]
                    label=""
                    if st_value!=1:
                        label=str(st_value)
                    G.add_edge(str(self.SpIdStrArray[j]), "r"+str(i), color="#4B0082",
                               label=label,title=label)
                
                for j in self.getIndArrayFromBt(self.ProdListBt[i]):
                    st_value=self.MpDf.iloc[j,i]
                    label=""
                    if st_value!=1:
                        label=str(st_value)
                    G.add_edge("r"+str(i), str(self.SpIdStrArray[j]), color="#4B0082",
                               label=label,title=label)
            
            elif disp_non_act:
                G.add_node("r"+str(i), color = "#E0E0E0", label="r"+str(i), shape="square", size=7,title="")
            
            
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
    
    def getCasualDecomGraphNx(self,decom):
        '''
        

        Parameters
        ----------
        decom : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
    
        # Creation of graph for causal decompostition
        Causal_G=nx.DiGraph()
       
        b_sp=bt(self.MpDf.shape[0])
        b_sp.setall(0)
    
        sp_EUF=b_sp.copy()
        for i in np.where(decom<0)[0]:
            sp_EUF[i]=1

        rp_EUF=self.getTriggerableRpBtFromSp(sp_EUF)
        # rp_X=self.getTriggerableRpBtFromSp(sp_set)
    
       
        # Assigning nodes to fragile cycles, overproduced and catalysts
        for i in np.unique(decom):
            
            if i>0:
                b_sp.setall(0)
                for k in np.where(decom==i)[0]:
                    b_sp[k]=1
                b_rp=self.getTriggerableRpBtFromSp(b_sp|sp_EUF)
                b_rp=b_rp&~rp_EUF
                Causal_G.add_node("D_"+str(int(i)), b_sp=b_sp.copy(), sp_id=self.SpIdStrArray[b_sp.search(1)],b_rp=b_rp,
                                  type="D")
            elif i==-1:
                
                j=1
                for k in np.where(decom==-1)[0]:
                    b_sp.setall(0)
                    b_sp[k]=1
                    
                    Causal_G.add_node("F_"+str(int(j)), b_sp=b_sp.copy(), sp_id=self.SpIdStrArray[b_sp.search(1)],
                                      b_rp=self.getTriggerableRpBtFromSp(b_sp),type="F")
                    j+=1
            
            elif i==-2:
                
                j=1
                for k in np.where(decom==-2)[0]:
                    b_sp.setall(0)
                    b_sp[k]=1
                    
                    Causal_G.add_node("E_"+str(int(j)), b_sp=b_sp.copy(), sp_id=self.SpIdStrArray[b_sp.search(1)],
                                      b_rp=self.getTriggerableRpBtFromSp(b_sp),type="E")
                    j+=1
                    
        
        sp_E=b_sp.copy()
        for i in np.where(decom==-2)[0]:
            sp_E[i]=1
        
        supp_sp=b_sp.copy()
        prod_sp=b_sp.copy() 
        print(Causal_G.nodes(data=True))
        # Assigning causality edges between nodes.
        for i in Causal_G.nodes(data=True):
            
            if i[1]['type']=="D":
                
                supp_sp.setall(0)
                prod_sp.setall(0)
                
                for j in i[1]['b_rp'].search(1):
                    supp_sp|=self.ReacListBt[j]
                    prod_sp|=self.ProdListBt[j]
                
                N_sp=supp_sp&~i[1]['b_sp']
                P_sp=prod_sp&~i[1]['b_sp']
                E_sp=sp_E&~i[1]['b_sp']
                
                P_nodes = list(map(lambda x: x[0],filter(lambda x: (x[1]['b_sp']&P_sp).any(), Causal_G.nodes(data=True))))
                N_nodes = list(map(lambda x: x[0],filter(lambda x: (x[1]['b_sp']&N_sp).any(), Causal_G.nodes(data=True))))
                
                for j in P_nodes:
                    if not Causal_G.has_edge(j, i[0]):
                        Causal_G.add_edge(j, i[0])
                
                for j in N_nodes:
                    if not Causal_G.has_edge(i[0],j):
                        Causal_G.add_edge(i[0],j)
                
                
            if i[1]['type']=="F":
                
                supp_sp.setall(0)
                prod_sp.setall(0)
                
                for j in i[1]['b_rp'].search(1):
                    supp_sp|=self.ReacListBt[j]
                    prod_sp|=self.ProdListBt[j]
                
                #N_sp=supp_sp
                P_sp=prod_sp
                # N_sp=supp_sp&~i[1]['b_sp']
                # P_sp=prod_sp&~i[1]['b_sp']
                
                P_nodes = list(map(lambda x: x[0],filter(lambda x: (x[1]['b_sp']&P_sp).any(), Causal_G.nodes(data=True))))
               # N_nodes = list(map(lambda x: x[0],filter(lambda x: (x[1]['b_sp']&N_sp).any(), Causal_G.nodes(data=True))))
                
                for j in P_nodes:
                    if not Causal_G.has_edge(j, i[0]):
                        Causal_G.add_edge(j, i[0])
                
                # for j in N_nodes:
                #     if not Causal_G.has_edge(i[0],j):
                #         Causal_G.add_edge(i[0],j)
                

                    
            if i[1]['type']=="E":
                
                supp_sp.setall(0)
                prod_sp.setall(0)
                
                for j in i[1]['b_rp'].search(1):
                    supp_sp|=self.ReacListBt[j]
                    prod_sp|=self.ProdListBt[j]
            
                N_sp=supp_sp
                P_sp=prod_sp
                # N_sp=supp_sp&~i[1]['b_sp']
                # P_sp=prod_sp&~i[1]['b_sp']
                
                P_nodes = list(map(lambda x: x[0],filter(lambda x: (x[1]['b_sp']&P_sp).any(), Causal_G.nodes(data=True))))
                N_nodes = list(map(lambda x: x[0],filter(lambda x: (x[1]['b_sp']&N_sp).any(), Causal_G.nodes(data=True))))
                
                for j in P_nodes:
                    if not Causal_G.has_edge(j, i[0]):
                        Causal_G.add_edge(j, i[0])
                
                for j in N_nodes:
                    if not Causal_G.has_edge(i[0],j):
                        Causal_G.add_edge(i[0],j)
                
        return Causal_G
                
    # # Function that generates the polyhedra and polytopes that the fragile and 
    # # overproduced cycles of a decomposition. It takes as an increment the 
    # # overproduced species opsp_set, the total species sp_set, and the present 
    # # reactions pr. It returns as objects the rays and vertices of the 
    # # overproduced set op_ph and lists them with the respective rays 
    # # and vertices of the fragile circuits fc_ph_list.   
    # def polyth(self,opsp_set,sp_set,pr):
    #     # Checks if input is a bitarray, if it is not, it make the 
    #     # transformation
    #     if not (isinstance(sp_set,bt)):
    #         sp=bt(self.MpDf.shape[0])
    #         sp.setall(0)
            
    #         for i in sp_set:
    #             if i in self.MpDf.index.values:
    #                 ind=self.MpDf.index.get_loc(i)
    #                 sp[ind]=1
    #     else:
    #         sp=sp_set
        
    #     if not (isinstance(opsp_set,bt)):
    #         opsp=bt(self.MpDf.shape[0])
    #         opsp.setall(0)
            
    #         for i in opsp_set:
    #             if i in self.MpDf.index.values:
    #                 ind=self.MpDf.index.get_loc(i)
    #                 opsp[ind]=1
    #     else:
    #         opsp=opsp_set
            
    #     if not (isinstance(pr,bt)):
    #         v=bt(self.MpDf.shape[0])
    #         v.setall(0)
            
    #         for i in pr:
    #             v[i]=1
    #     else:
    #         v=pr
        
        
        
    #     Ns=sp.count() # number of species (rows)
    #     nsp=list(set(range(len(sp)))-set(self.getIndArrayFromBt(sp))) #species no present
        
        
    #     rc =[] # creating variable of available reaction 
    #     for j in self.getIndArrayFromBt(v):
    #         if (all(self.MpDf.iloc[nsp,j]==0) and all(self.MrDf.iloc[nsp,j]==0)):
    #             rc.append(j) # selecting reactions that can be trigger with available species
        
    #     E = np.identity(len(v))
    #     f = np.zeros(len(v))
        
    #     proj =(E,f) 
       
    #     # stoichiometric matrix of rn with prepended destruction reactions
    #     S=self.MpDf-self.MrDf
    #     S=S.to_numpy()  
        
    #     S=np.concatenate((S,np.identity(len(v))),axis=0)
        
    #     p=np.zeros(Ns+len(v))
    #     p[self.getIndArrayFromBt(opsp)]=1
        
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
            
            
    #     dcom=self.getDcomArray(opsp,sp,v)
    #     print("dcom :",dcom)
    #     fc_ph_list=[]
        
    #     for i in np.unique(dcom[dcom>0]):
            
    #         print("i: ",i)
    #         sp_ind=np.where(dcom==i)[0]
    #         # print("sp_ind: ",sp_ind)
            
    #         # print("mr condition: ",self.MrDf.iloc[sp_ind,:].sum(axis=0)>0)
    #         # print("mp condition: ",self.MpDf.iloc[sp_ind,:].sum(axis=0)>0)
            
    #         sub_rc=set(np.where(self.MrDf.iloc[sp_ind,:].sum(axis=0)>0)[0])
    #         sub_rc=set(sub_rc) | set(np.where(self.MpDf.iloc[sp_ind,:].sum(axis=0))[0])
    #         sub_rc=np.array(list(sub_rc))
   
    #         sub_sp=set(np.where(self.MrDf.iloc[:,sub_rc].sum(axis=1)>0)[0])
    #         sub_sp=set(sub_sp) | set(np.where(self.MpDf.iloc[:,sub_rc].sum(axis=1)>0)[0])
    #         # sub_sp=np.array(list(sub_sp))
            
    #         print("sub_sp: ",sub_sp)
            
    #         nsp=list(set(range(len(sp)))-sub_sp) #species no present
    #         sub_rc =[] # creating variable of available reaction 
    #         for j in self.getIndArrayFromBt(v):
    #             if (all(self.MpDf.iloc[nsp,j]==0) and all(self.MrDf.iloc[nsp,j]==0)):
    #                 sub_rc.append(j) # selecting reactions that can be trigger with available species
           
    #         print("sub_rc: ",sub_rc)
            
    #         op_ind=np.array(set(self.getIndArrayFromBt(opsp)) & set(sp_ind))
    #         cy_ind=np.array(set(dcom[dcom==-2]) & set(sp_ind))
            
    #         E = np.identity(len(v))
    #         f = np.zeros(len(v))
            
    #         proj =(E,f)
        
    #         S=self.MpDf-self.MrDf
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
            
            
                

