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
from scipy.optimize import linprog



file="../RN_software/rn_test.txt"


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
                        st=re.search("[0-9]*",rn[2][i][j]).group()
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
            out.mr=mr
            out.mp=mp
            out.reac=reac
            out.prod=prod
            
            
            return out
        except:
            print("error reading smbl file")
    # Initialization of the network from a smbl file, "file" corresponds 
    # to the path of the smbl file. 
    @classmethod
    def from_sbml(cls,file):
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
                
                if not (i.find('reversible') is None):
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
            out.mr=mr
            out.mp=mp
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
            sp=bt(self.mp.shape[1])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.columns.values:
                    ind=self.mp.columns.get_loc(i)
                    sp[ind]=1
        else:
            sp=copy.copy(sp_set)
        # creating a vector of reaction that can be trigered
        n_reac = np.array(range(self.mp.shape[0]))
        i=0
        flag=False
        
        # Generates the closure until no reaction can be trigered
        while(len(n_reac)>0):
            if (sp & self.reac[n_reac[i]]).count() == self.reac[n_reac[i]].count():
                
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
            sp_set=np.array(range(self.mp.shape[1]))
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
            sp=bt(self.mp.shape[1])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.columns.values:
                    ind=self.mp.columns.get_loc(i)
                    sp[ind]=1
                      
        # Generates the closure only for the first pass
        for i in range(self.mp.shape[0]):
            if (sp & self.reac[i]).count() == self.reac[i].count():
                sp=sp|self.prod[i]
        
        # returns bitarray or species vector
        if bt_type : 
            return sp
        else:
            sp_set=np.array(range(self.mp.shape[1]))
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
            sp=bt(self.mp.shape[1])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.columns.values:
                    ind=self.mp.columns.get_loc(i)
                    sp[ind]=1
                      
        # init of reactant species and product species array
        r = bt(self.mp.shape[1])
        r.setall(0)
        
        p = bt(self.mp.shape[1])
        p.setall(0)
        
        # generates of reactant species and product species array
        for i in range(self.mp.shape[0]):
             if (sp & self.reac[i]).count() == self.reac[i].count():
                r|=self.reac[i]
                p|=self.prod[i]
        
        # verifies if produces species are in the sp set and if ssm condition 
        # is statisfy
        if (sp & p).count() == p.count():
            return ((r & p).count() == r.count())
        else:
            return False

    # Function that confirms if a set is stoichimetriclly semi-self-mantained, 
    # input is "sp_set" that can be an bitarray or an species array. 
    # returns a true or false depending if the property is achived
    def is_s_ssm(self,sp_set):
        
        # Checks if input is or not bita array, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[1])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.columns.values:
                    ind=self.mp.columns.get_loc(i)
                    sp[ind]=1
                      
        # init of reactant species and product species array
        r = bt(self.mp.shape[1])
        r.setall(0)
        
        p = bt(self.mp.shape[1])
        p.setall(0)
        
        # generates of reactant species and product species array
        for i in range(self.mp.shape[0]):
             if (sp & self.reac[i]).count() == self.reac[i].count():
                for j in range(self.mp.shape[1]):
                    d=self.mp.iloc[i][j]-self.mr.iloc[i][j]
                    if d<0:
                        r[j]=1
                    elif d>0:
                        p[j]=1
        
        # verifies if produces species are in the sp set and if ssm condition 
        # is statisfy
        if (sp & p).count() == p.count():
            return ((r & p).count() == r.count())
        else:
            return False

    def is_sm(self,sp_set):
        
        # Checks if input is or not bita array, if it's no, it make the 
        # transmation
        if not (isinstance(sp_set,bt)):
            sp=bt(self.mp.shape[1])
            sp.setall(0)
            
            for i in sp_set:
                if i in self.mp.columns.values:
                    ind=self.mp.columns.get_loc(i)
                    sp[ind]=1
                    
        # init of reactant species and product species array
        r = bt(self.mp.shape[1])
        r.setall(0)
        
        p = bt(self.mp.shape[1])
        p.setall(0)
        
        # generates of reactant species and product species array
        for i in range(self.mp.shape[0]):
             if (sp & self.reac[i]).count() == self.reac[i].count():
                r|=self.reac[i]
                p|=self.prod[i]
        
        # adding product species to the species vector
        sp=r|p
        
        # init of reaction indexes form trigable reactions
        r_ind=[]
        
        # generates of reaction indexes form trigable reactions
        for i in range(self.mp.shape[0]):
             if (sp & self.reac[i]).count() == self.reac[i].count():
                 r_ind.append(i)
        
        # init of  indexes form trigable reactions
        sp_ind=[]
        for i in range(len(sp)):
            if sp[i]==1:
                sp_ind.append(i)
        
        # relative Soichiometric matrix to present species and reaction
        S=np.array(self.mr-self.mp)
        S=np.transpose(S[np.ix_(r_ind,sp_ind)])
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
    
    # Function of species that returns the partiton that contains the species
    def sp2p(self, sp):
        p=bt(len(self.sp_b))
        p.setall(0)
        for i in range(len(self.sp_b)):
            if (sp & self.sp_b[i]).count() == self.sp_b[i].count():
                p[i]=1
        
        return p
    
    # Subset funcntion for bitarrays
    def is_subset(a, b):
        return (a & b).count() == a.count()
      
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
        x_r_p=-np.ones(self.mp.shape[0])
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
            
            for j in range(i,self.mp.shape[0]):
                if c_reac[i]==c_reac[j]:
                    x_r_p[j]=xeqc
                    
            xeqc+=1
        
        self.sp_b=sp_b
        self.x_r_p=x_r_p
        
        # assingnig species contained each partition
        sp_p=[]
        
        b_sp=bt(self.mp.shape[1])
        b_sp.setall(0)
        
        for i in range(xeqc):
            sp_p.append(b_sp.copy())
            
        for i in range(self.mp.shape[0]):
            sp_p[x_r_p[i]]|= self.reac[i]|self.prod[i] 
        
        
        self.sp_p=sp_p
        
        
        # assingnig reaction contained each partition
        r_p=[]
        b_r=bt(self.mp.shape[0])
        b_r.setall(0)
        
        for i in range(xeqc):
            r_p.append(b_r)
        for i in range(self.mp.shape[0]):
            r_p[x_r_p[i]][i]=1
            
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
                    psp_b[i]|=self.prod[k];
                    sp_sp_b[i]|=bt(self.mr.iloc[k,:] < self.mp.iloc[k,:])
                    sn_sp_b[i]|=bt(self.mr.iloc[k,:] > self.mp.iloc[k,:])
            
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
                if (not j==i) and p_b[j][j]==0 and (rsp_b[i] & psp_b[j]).any():
                    dyn_conn[i][j]=1
                if (not j==i) and p_b[j][j]==0 and (psp_b[i] & rsp_b[j]).any():
                    dyn_conn[j][i]=1
                # Hasse condition (all conected)
                if (not j==i) and p_b[i][j]==0:
                    conn[i][j]=1
                    conn[j][i]=1
        
        self.conn=conn
        self.dyn_conn=dyn_conn
    
    
    #for a given bitset of containend basic set, return the connected basics 
    # sets to it
    def conn_b(self,s):
        c=bt(len(self.sp_b))
        c.setall(0)
        
        for i in self.bt_ind(s):
            c|=self.conn[i]
            
        c-=s
        return c
        
    #for a given bitset of containend basic set, return the dynamically 
    # connected basics sets to it
    def dyn_conn_b(self,s):
        c=bt(len(self.sp_b))
        c.setall(0)
        
        for i in self.bt_ind(s):
            c|=self.dyn_conn[i]
            
        c-=s
        return c
    
    #for a given bitset of containend basic set, return the basic sets that 
    # can contribute so the semi-self-maintainace can be reached
    def contrib_b(self,s):
        
        # negative sotichometry species
        
        n=bt(len(self.mp.shape[1]))
        n.setall(0)
        p=n.copy()
        pp=n.copy()
        
        for i in self.bt_ind(s):
            p|=self.sp_sp_b[i]
            n|=self.sn_sp_b[i]
                       
        # P correspond to bitarray whit positive stoichometry
        n=-p
        
        c=s.copy()
        if n.empty():
           return c.setall(0)
       
        # posible sets that can contribute to be semi-self-mantained
        c.flip()
        
        # eliminatig basics form c that no will contrubute
        for i in self.bt_ind(c):
            if not (self.sp_sp_b[i] & n).any():
                c[i]=0
            else:
                pp|=self.sp_sp_b[i]
        
        # pp are posible produces species
        # if they don fullfil consumed species then ther is no contribution
        # to be a semi-self-mantained
        if not (n-pp).empty():
            return c.setall(0)
        else:
            return c
        
        
                
    ´