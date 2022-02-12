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
    # Example rn = RN.("rn_text.txt")
    
    @classmethod
    def form_txt(cls,file):
        try:
            # list fitering
            rn = pd.read_csv(file,header=None)
            rn[0] = rn[0].str.replace("[ \t]","",regex=True)        
            rn[0] = rn[0].str.replace("^.*:","",regex=True)        
            rn[0] = rn[0].str.replace("[;#].*$","",regex=True)
            rn[1] = rn[0].str.contains("=>")         
            rn[2] = rn[0].str.split("->|=>",expand=True)[0]
            rn[3] = rn[0].str.split("->|=>",expand=True)[1]
            rn[2] = rn[2].str.split("+",expand=False)
            rn[3] = rn[3].str.split("+",expand=False)
            
            print(rn)
            # obtaining list of reactions separated as reactive part
            er = list()
            ep = list()
            
            for i in range(len(rn[0])):
                # print(i)
                psp = list()
                pst = list()
                ssp = list()
                sst = list()
               
                
                if len(rn[2][i])!=0:
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
                
                if len(rn[3][i])!=0:
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
                sp=set.union(sp,set(er[i][0]),set(ep[i][0]))
            
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
            
            
            return out
        except:
            print("error reading txt file")
    # Initialization of the network from a smbl file, "file" corresponds 
    # to the path of the smbl file. 
    @classmethod
    def from_sbml(cls,file,modifiers=True):
        try:
            with open(file,'r') as file:
                bs_sbml = bs(file.read(), "xml")
        
            er = list()
            ep = list()

            b_i=bs_sbml.findAll(boundaryCondition="true")
        
            for i in b_i:

                # print(i['id'])
                er.append([[[]],[[]]])
                ep.append([[i['id']],[1]])
            
            # obtaining list of reactions separated as reactive part 
            # and product part
            bs_r=bs_sbml.select('reaction[id]')            

    
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
                    # print(i['id'])
                    reac = reac.select('speciesReference')
                    for j in reac:

                        ssp.append(j['species'])
                        if j.has_attr('stoichiometry'):
                            sst.append(float(j['stoichiometry'])) 
                        else:
                            sst.append(1.0)
                else:
                        ssp.append([])
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
                        psp.append([])
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
            
            return out
        except:
            print("error reading smbl file")
            
    
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
                # print(sp)
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
        
        # print(sp)
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
        
        bounds=[] #creation of bounds vectror, solution sould be greater that 1 
        # for each component
        for i in range(len(r_ind)):
            bounds.append((1,None))

        # lineal programing calculation of existance of a solution to be 
        # self-mantainend
        
        res = linprog(c, A_ub=S, b_ub=b, bounds=bounds)
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
        
        # creating equivalence clases for each reaction that generate the same
        # closure and so each basic can be created
        for i in range(len(x_r_p)):
            # reaction already assigned to equivalnce class
            if x_r_p[i]>=0:
               continue
           
            x_r_p[i]=xeqc
            sp_b.append(c_reac[i])
            
            for j in range(i,self.mp.shape[1]):
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
        
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.p_b)):
            G.add_node(fbt(self.p_b[i]),level=self.p_b[i].count(),
                       sp=self.sp[self.bt_ind(self.sp_b[i])],
                       is_ssm=self.is_ssm(self.sp_b[i]))
            if self.is_ssm(self.sp_b[i]):
                ssms.append(self.sp[self.bt_ind(self.sp_b[i])])
        
        # Generation of the multigraph of the synergic structure by set level (contained basics)
        for i in range(len(self.sp)):
             
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
             
             # Generating closrues whit connected basics sets for each set in level i
            for j in nodes:
                 for k in self.bt_ind(self.conn_b(bt(j))):
                     
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
        self.ssms=ssms
                    
                    
       
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
        
        # The nodes corresponding to the basic sets are generated.
        for i in range(len(self.p_b)):
            G.add_node(fbt(self.p_b[i]),level=self.p_b[i].count(),
                       sp=self.sp[self.bt_ind(self.sp_b[i])],
                       is_ssm=self.is_ssm(self.sp_b[i]))
            if self.is_ssm(self.sp_b[i]):
                ssms.append(self.sp[self.bt_ind(self.sp_b[i])])
        
        # Generation of the multigraph of the synergic structure by set level (contained basics)
        for i in range(len(self.sp)):
             
            # Closed set (nodes) at level i
            nodes = [x for x,y in G.nodes(data=True) if y['level']==i+1]
             
             # Generating closrues whit connected basics sets for each set in level i
            for j in nodes:
                
                # if node is semi-self-maintained (ssm), the it explore other possible 
                # combinations to serch for ssm sets, if not search for basic that can 
                # contibubte to be ssm
                
                if self.syn_str.nodes[j]["is_ssm"]:
                     conn=self.bt_ind(self.conn_b(bt(j)))
                else:
                     contib=self.contrib_b(bt(j))
                     if not contib.any():
                         continue
                     conn=self.bt_ind(contib)
                 
                for k in conn:
                     
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
        self.ssms=ssms
    
                
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
                    
       
    # verifies which species of a reaction network (not necessarily an organization) 
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
        
        # print("S: ",np.array(S))

        
        
        # overproducible status for species (initially False, until proven to be True).
        o=[]
        for i in range(len(sp)):
            o.append(False)
        
        for p in range(len(sp)):
            if not (o[p]):  # already known overproducible species are skipped 
                
                S[p][p]=1 # creation instead of destruction for p
                f[p]=1 # production of p = 1 instead of 0, if it is possible with creation rate 0, then it is overproducible
                cost[p]=1  # cost for creation of p = 1 instead of 0
            
                # lineal programing calculation of existance of a solution to be 
                # self-mantainend
                
                print("S: ",np.array(S))
                print("p: ",p)
                print("f: ",f)
                print("cots: ",cost)
               
                res = linprog(cost, A_eq=S, b_eq=f,method='highs')
                # The unknown is the process vector v. The equations are S v = f with the inequality constraint v>=0 (rates can be 0).
                # Only the creation reaction for p is penalized in Cost (the rate should be 0 if p is overproducible).
                print("x: ",res.x)
                print("\n")
                
                if(res.x[p]==0):
                
                    o[p] = True # no need of creation implies p is overproducible
                    for i in range(Ns):
                        if res.x[i]>0:
                            o[i]=True
            
                S[p][p] = -1 
                f[p] = 0
                cost[p] = 0  # original destruction reaction and zero values for next iteration
            
        
        print(o)
        opsp=bt(len(sp))
        opsp.setall(0)
        for i in range(len(o)):
            if o[i]:
                opsp[self.bt_ind(sp)[i]]=1
                
        return(opsp)
    
    