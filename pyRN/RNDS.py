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
import pypoman as ph
from itertools import chain, combinations
from typing import List, Any, Iterable


class RNDS(RNIRG):
    
    # Verifies which species of a reaction network (not necessarily an organization) 
    # are overproducible. Inputs are species sp_set and process vector pr, 
    # returns a list of overproducible species
    def over_prod(self,sp_set,pr):     
        # Checks if input is bitarray, if it is not, it make the 
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
        # Checks if input is bitarray, if it is not, it make the 
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
        
        # Checks if input is bitarray, if it is not, it make the 
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
                    
                                       
    # Function that verifiys if a decomposition input is semi-self-
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
        # Checks if input is a bitarray, if it is not, it make the 
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
        # Checks if input is a bitarray, if it is not, it make the 
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
        # Checks if input is a bitarray, if it is not, it make the 
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
    #     # Checks if input is a bitarray, if it is not, it make the 
    #     # transformation
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
            
            
                

