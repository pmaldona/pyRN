#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 10:48:31 2022

@author: pmaldona

Reaction Network Simulator and Random Walk Class
"""
from .RNIRG import RNIRG
import numpy as np
import roadrunner as re
import pandas as pd
import time


class RNSRW(RNIRG):
    # Funtions that create a mass action dinamics telurrium model (CRNS.model) of 
    # reaction network. It recives as input a vector of the initial concentration 
    # of species i_sp, the reactive constant vector rt and the concentration 
    # thearshold where a species is condidered not present. if the variables 
    # i_sp and rt are not given, they are randomly initialized by a uniform 
    # distribution between 0 and 1
    def set_model_ma(self, i_sp=None ,rt=None ,cutoff=0.1):
        
        # Creating the model
        # self.model=te.loada("")
        try:
            self.model
            self.model.clearModel()
        except:
            self.model=re.RoadRunner()
        
        self.model.addCompartment("C", 1,True)
        # Creating the random initial concetration if it's out of condition
        if i_sp is None: 
            i_sp=np.random.rand(self.mp.shape[0])
        elif(len(i_sp) !=self.mp.shape[0]):
            i_sp=np.random.rand(self.mp.shape[0])
        
        # Creating the random initial rates if it's out of condition
        if rt is None:  
            rt=np.random.rand(self.mp.shape[1])
        elif (len(rt) !=self.mp.shape[1]):
           rt=np.random.rand(self.mp.shape[1])
           
        # Adding species to the model and treshold concetration
        for i in range(len(self.sp)):
            self.model.addSpecies(self.sp[i], compartment="C", initConcentration=i_sp[i],forceRegenerate=True)
        # Adding reactions and reactions rate constants
        for i in range(self.mp.shape[1]):
            self.model.addParameter("k"+str(i),rt[i],True)
            
            prod=[]
            for j in np.where(self.mp.iloc[:,i])[0]:
                if self.mp.iloc[j,i].astype(int) > 1:
                    prod.append(str(self.mp.iloc[j,i].astype(int))+self.mp.index[j])
                else:
                    prod.append(self.mp.index[j])
            reac=[]
            rate="k"+str(i)
            for j in np.where(self.mr.iloc[:,i])[0]:
                if self.mr.iloc[j,i].astype(int) > 1:
                    reac.append(str(self.mr.iloc[j,i].astype(int))+self.mr.index[j])
                else:
                    reac.append(self.mr.index[j])
                    
                if self.mr.iloc[j,i]>1:
                    for k in range(1,self.mr.iloc[j,i].astype(int)+1):
                        rate+=" * "+self.mr.index[j]
                else:
                    rate+=" * "+self.mr.index[j]
   
            self.model.addReaction("r"+str(i), reac, prod, rate)
        
        # Adding events so if conentrations of species are under cutoff, 
        # will have zero concentration 
        # for i in range(len(self.sp)):
        #     self.model.addEvent("e"+str(i),True,self.sp[i]+"<"+str(cutoff),True)
        #     self.model.addEventAssignment("e"+str(i),self.sp[i],"0",True)
        #     self.model.addTrigger("e"+str(i),self.sp[i]+"<"+str(cutoff),forceRegenerate=True)
        
        self.model.regenerateModel()
        self.sbml=False
        
    # Function that load dynamical model directly form the sbml file     
    def set_sbml_model(self):
        if not self.sbml:
            print("network do not correspond to an sbml file")
            return
        # self.model = te.loadSBMLModel(self.fname)
        self.model = re.RoadRunner(self.fname)
        
    
    # Function that simulates the dynamics of the proposed model. It takes as 
    # input the initial time ti and the final time tf, and the number of steps 
    # to simulate by means of the variable steps, and a cutoff concentration  
    # that considers if a species is present or not. The function generates the 
    # dataframes RNIRG.con with RNIRG.rate, which correspond to concentrations and 
    # processes of the simulation respectively. It also generates tre more variables, 
    # RNIRG.abst, RNIRG.a_sp and RNIRG.a_r. Those corresponde to the abstraction, 
    # active species and active reaction. This dataframes will append another 
    # line if any of the species change their concentration below or over the cutoff  
    def run_model(self,ti=0,tf=50,steps=100,cutoff=0.1): 
        
        
        # Generation of kinetic contstant vector
        k=[]
        rate_name=[]
        for i in range(self.mr.shape[1]):
            k.append(self.model.getGlobalParameterByName('k'+str(i)))
            rate_name.append('r'+str(i+1))
        
        k=np.array(k)
        
        # Considering if the simulation has benn aready run, and generation of
        # output variables
        

        t_step=(tf-ti)/steps
        c_st=self.model.getFloatingSpeciesConcentrationsNamedArray()[0].copy()
        c_st[c_st<cutoff]=0
        con=[[ti]+c_st.tolist()]
        rate=[[ti]+self.model.getReactionRates().tolist()]
        abst=[[ti]+self.abstract_sp(self.model.getFloatingSpeciesConcentrationsNamedArray()[0],k).tolist()]
        a_sp=[[ti]+self.active_sp(self.model.getFloatingSpeciesConcentrationsNamedArray()[0]).tolist()]
        a_r=[[ti]+self.active_reac(self.model.getFloatingSpeciesConcentrationsNamedArray()[0]).tolist()]
        ac_reac=self.active_reac(self.model.getFloatingSpeciesConcentrationsNamedArray()[0])
        # Running the simulation
        for i in range(steps):
            t=t_step*i+ti
            self.model.oneStep(t,t_step)
            n_st=self.model.getFloatingSpeciesConcentrationsNamedArray()[0].copy()
            n_st[n_st<cutoff]=0
            con.append([t+t_step]+n_st.tolist())
            reac=self.model.getReactionRates()
            reac[ac_reac==False]=0
            rate.append([t+t_step]+reac.tolist())
            # rate.append([t+t_step]+self.model.getReactionRates().tolist())
            # appending data if the active species changes
            if not all((c_st<cutoff)==(n_st<cutoff)):
                abst.append([t+t_step]+self.abstract_sp(self.model.getFloatingSpeciesConcentrationsNamedArray()[0],k).tolist())
                a_sp.append([t+t_step]+self.active_sp(self.model.getFloatingSpeciesConcentrationsNamedArray()[0]).tolist())
                ac_reac=self.active_reac(self.model.getFloatingSpeciesConcentrationsNamedArray()[0])
                a_r.append([t+t_step]+ac_reac.tolist())
            c_st=n_st.copy()
            c_st[c_st<cutoff]=0
        
        
        # Generating output
        con=pd.DataFrame(con,columns=['time']+self.sp.tolist())
        con=con.set_index('time')
        
        abst=pd.DataFrame(abst,columns=['time']+self.sp.tolist())
        abst=abst.set_index('time')
        
        a_sp=pd.DataFrame(a_sp,columns=['time']+self.sp.tolist())
        a_sp=a_sp.set_index('time')
        
        
        rate=pd.DataFrame(rate,columns=['time']+rate_name)
        rate=rate.set_index('time')
        
        a_r=pd.DataFrame(a_r,columns=['time']+rate_name)
        a_r=a_r.set_index('time')
        
        
        try:
            self.con=pd.concat([self.con,con])
            self.rate=pd.concat([self.rate,rate])
            self.abst=pd.concat([self.abst,abst])
            self.a_sp=pd.concat([self.a_sp,a_sp])
            self.a_r=pd.concat([self.a_r,a_r])
            
        except:
            self.con=con
            self.rate=rate
            self.abst=abst
            self.a_sp=a_sp
            self.a_r=a_r

    
    
    # Function that modify parameters of a given model. It recives as input a 
    # vector of the initial concentration of species i_sp, the reactive 
    # constant vector rt and the concentration 
    def param_change_ma(self, i_sp=None ,rt=None,init=False):
        
        if self.sbml:
            print("The initialized model corresponds to an sbml model, therefore the perturbation function cannot be used.")
            return
        
        # Creating the random initial concetration if it's out of condition
        if i_sp is None: 
            i_sp=np.random.rand(self.mp.shape[0])
        elif(len(i_sp) !=self.mp.shape[0]):
            i_sp=np.random.rand(self.mp.shape[0])
        
        # Creating the random initial concetration if it's out of condition
        if rt is None:  
            rt=np.random.rand(self.mp.shape[1])
        elif (len(rt) !=self.mp.shape[1]):
            rt=np.random.rand(self.mp.shape[1])
        
        if init:
            del(self.con)
            del(self.rate)
            del(self.abst)
            del(self.a_r)
            del(self.a_sp)

            for i in range(self.mp.shape[1]):
                self.model.setGlobalParameterByName("k"+str(i),rt[i])
            
            for i in range(self.mp.shape[0]):
                self.model.setInitConcentration(self.model.getIds()[i], i_sp[i])
                
        else:
            
            for i in range(self.mp.shape[1]):
                self.model.setGlobalParameterByName("k"+str(i),rt[i])
            
            for i in range(self.mp.shape[0]):
                self.model.setValue(self.model.getIds()[i], i_sp[i])
    
    # randomize components of a state or flow vector using a log normal distribution
    # v vector, mask of affected components (defaults to active ones), mu and sigma parameters of the distribution
    def pert_randomize(self,v,mask=None,mu=0,sigma=.5):
        if mask is None:
            mask=(v!=0)
        mask = np.repeat(True,len(v)) & mask
        k = np.where(mask)[0]
        if (len(k)==0): 
            return(v)
        v[k] = np.exp(np.random.normal(mu,sigma,len(k)))
        return v
    
    # activation of components of v (changing values from 0 to 1 or from >0 to 0), mask (defaults to all components) 
    # n total active components required, up to p active components to be preserved (defaults to all active components)
    def pert_activation(self,v,mask=None,n=None,p=None):
        if mask is None:
            mask=np.repeat(True,len(v))
        if n is None:
            n = np.sum(mask)
        if p is None:
            p = np.sum(v[mask]>0)
        
        w = v[mask]
        l = len(w)
        a = np.sum(w>0)
        n = np.max([np.min([n,l]),0])
        p = np.max([np.min([p,a]),0])
        if n==0: 
            w = 0*w
        elif n<=p:
            if a==1: 
                k = np.where(w>0)[0] 
            else: 
                k = np.random.choice(np.where(w>0)[0],n)+1
            w[-k] = 0
        elif p==0:
            if l==1: 
                k = 1 
            else: 
                k = np.random.choice(l,n)+1
            w[-k] = 0
            w[k[w[k-1]==0]-1] = 1
        else :  # n>p and p>0
            i = np.where(w>0)[0]
            if len(i)>1: 
                i = np.random.choice(i,p)+1  # the ones we keep
            j = np.array(range(l))[-i]  # the remaining components
            k = np.random.choice(len(j),n-p)+1
            w[j[-k]] = 0
            k <- j[k-1]-1
            w[k[w[k-1]==0]-1] = 1
  
        v[mask] = w
        return v

    # perturb a vector by adding or substracting components and randomizing the result
    # v vector, d component delta, at least nmin active components, sigma random dispersion
    def pert_delta(self,v,d=1,nmin=5,sigma=.5):
        n = np.max([np.sum(v>0)+d,nmin])
        v = self.pert_activation(v,n=n)
        v = self.pert_randomize(v,sigma=sigma)
        
        return v


    # creates a structure to study evolutive paths of a reaction network rn
    # result: the reaction network rn, a random walk list rw
    # each rw[[i]] contains matrices f, s, p, c, a and r were each column stores
    # the starting state (f,s), the perturbed state p, its convergence (c), its abstraction (a) and used species (u)
    # for each perturbation and simulation step applied
    def rw_start(self):
      self.rw=[]
      self.rw.append(dict(f=None,s=None,p=None,c=None,a=None,u=None))
    
    
    # returns the boolean abstract state (closure) for a given reaction network rn
    # species concentrations vector s and flow constants vector f
    def abstract_sp(self,s,f,cutoff=.1):
        k = self.bt_ind(self.closure(sp_set=self.sp[np.where(s>=cutoff)[0]],r_set=np.where(f>0)[0],bt_type=True))
        v = np.repeat(False,len(s))
        v[k] = True
        return v
    
    # returns the boolean active reaction
    # species concentrations vector s
    def active_reac(self,s,cutoff=.1):
        sp = np.where(s>cutoff)[0]        
        def ac_reac(reac):
            return all(np.isin(np.where(reac)[0],sp))
        
        k=np.where(np.apply_along_axis(ac_reac, 0, self.mr))[0]
        v = np.repeat(False,self.mr.shape[1])
        v[k] = True
        return v


    # returns the boolean active species 
    # species concentrations vector s
    def active_sp(self,s,cutoff=.1):
        sp = np.where(s>cutoff)[0]
        def ac_reac(reac):
            return all(np.isin(np.where(reac)[0],sp))
        
        v = np.repeat(False,len(s))
        reac = np.where(np.apply_along_axis(ac_reac, 0, self.mr))[0]
        if len(reac)>0:
            sp = np.where((self.mr.iloc[:,reac].sum(axis=1)+self.mp.iloc[:,reac].sum(axis=1))>0)[0]
            v[sp] = True
            
        return v
    
    
    # working script to generate a random reaction network and perturbation/simulation random walks
    # if e is provided new random walks or steps are added to e
    # if rn is provided that reaction network is used instead of generating a random one
    # the random walks to be created or completed are defined by range w, by default 1:10
    # by default the number of perturbation and simulation steps l is set to 10
    # cutoff is the concentration threshold for a species to be reactive and n is the number of simulation steps 
    # the result is available in global variable e (the value is also returned by the function)
    def scr_gen_and_pert(self,w=range(10),l=10,cutoff=.1,n=5000):
                
        try:
            self.rw
        except:        
            self.rw_start() # the starting structure to store the random walks to be generated
                
        for i in w:  # for each random walk
            if i>len(self.rw)-1: # this is a new random walk
                    self.rw.append(dict(f=None,s=None,p=None,c=None,a=None,u=None,t=None)) # matrices are created to store the steps in columns
            if self.rw[i]['f'] is None: # this is a void random walk (0 steps)
                s = np.zeros(len(self.sp)) # the current state is zeroed
                f = self.pert_randomize(np.ones(self.mr.shape[1])) # the flow vector is randomized (each random walk has a different f)
                if i==0:
                    self.set_model_ma(i_sp=s ,rt=f ,cutoff=cutoff) # initialization of the model
                else:
                    self.param_change_ma(i_sp=s ,rt=f,init=True)
                    
            else: # this is a random walk with a number of steps already accumulated
              j = self.rw[i]['c'].shape[1]  # the number of steps up to now
              s = self.rw[i]['c'][:,j] # the exploration continues from the last convergence state in the random walk
              f = self.rw[i]['f'][:,j] # the last flow vector is conserved
            
            for k in range(l):  # for each j step in random walk i
                print("walk: "+str(i+1)+", step: "+str(k+1))
                # the current state is stored in the random walk
                
                if k==0:
                    s = self.pert_delta(np.zeros(len(self.sp)),d=1,nmin=1,sigma=.5)
                    self.rw[i]['s']=pd.DataFrame(self.model.getFloatingSpeciesConcentrationsNamedArray(),columns=self.sp).T
                    self.param_change_ma(s)
                else:
                    self.rw[i]['s'] = pd.concat([self.rw[i]['s'],self.con.iloc[-1]],axis=1)
                    s = self.pert_delta(np.array(self.con.iloc[-1]),d=1,nmin=1,sigma=.5) # a delta perturbation is applied to current state
                    self.param_change_ma(s)
                # end perturbations, start simulation:
                start = time.time() 
                self.run_model(n*k,n*(k+1),n,cutoff) # the perturbed state is simulated reaching a convergence state
                end = time.time()
                st=end-start


                # end simulation
                if k==0:
                    self.rw[i]['f'] = self.rate.iloc[-1].copy() # the flow vector is stored in the random walk
                    self.rw[i]['p'] = pd.DataFrame(s,index=self.sp)  # the perturbed current state is stored in the random walk
            # if (is.null(cs)) cs <- s # if something went wrong with the simulation, we just keep the initial current estate
                    self.rw[i]['c'] = self.con.iloc[-1].copy() # the convergent state is stored in the random walk
                    self.rw[i]['a'] = self.abst.iloc[-1].copy() # the abstraction is stored  
                    self.rw[i]['u'] = self.a_r.iloc[-1].copy() # a second abstraction is stored (used species)
                    self.rw[i]['t'] = [st] # the time elapsed in the dynamic simulation is stored in the random walk
                else:
                    self.rw[i]['f'] = pd.concat([self.rw[i]['f'],self.rate.iloc[-1]],axis=1) # the flow vector is stored in the random walk
                    self.rw[i]['p'] = pd.concat([self.rw[i]['p'],pd.DataFrame(s,index=self.sp)],axis=1)  # the perturbed current state is stored in the random walk
              # if (is.null(cs)) cs <- s # if something went wrong with the simulation, we just keep the initial current estate
                    self.rw[i]['c'] = pd.concat([self.rw[i]['c'],self.con.iloc[-1]],axis=1) # the convergent state is stored in the random walk
                    self.rw[i]['a'] = pd.concat([self.rw[i]['a'],self.abst.iloc[-1]],axis=1) # the abstraction is stored  
                    self.rw[i]['u'] = pd.concat([self.rw[i]['u'],self.a_r.iloc[-1]],axis=1) # a second abstraction is stored (used species)
                    self.rw[i]['t'].append(st) # the time elapsed in the dynamic simulation is stored in the random walk
        
            self.rw[i]['f'].columns=range(self.rw[i]['f'].shape[1])
            self.rw[i]['s'].columns=range(self.rw[i]['s'].shape[1])
            self.rw[i]['p'].columns=range(self.rw[i]['p'].shape[1])
            self.rw[i]['c'].columns=range(self.rw[i]['c'].shape[1])
            self.rw[i]['a'].columns=range(self.rw[i]['a'].shape[1])
            self.rw[i]['u'].columns=range(self.rw[i]['u'].shape[1])