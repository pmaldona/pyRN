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
import json
import copy

class RNSRW(RNIRG):
    # Funtions that create a mass action dinamics telurrium model (CRNS.model) of 
    # reaction network. It recives as input a vector of the initial concentration 
    # of species i_sp, the reactive constant vector rt and the concentration 
    # thearshold where a species is condidered not present. if the variables 
    # i_sp and rt are not given, they are randomly initialized by a uniform 
    # distribution between 0 and 1
    def setMakModel(self, i_sp=None ,rt=None ,cutoff=0.1):
        
        self.CutoffFloat=cutoff
        # Creating the model
        # self.model=te.loada("")
        try:
            self.model
            self.model.clearModel()
            try:
                del(self.SpConDf)
                del(self.RpRateDf)
            except:
                self.model
        except:
            self.model=re.RoadRunner()
        # self.model.setIntegrator('rk45') # set integrator first
        # self.model.integrator.epsilon = 1e-10
        # self.model.getIntegrator().setValue("absolute_tolerance",[1E-18,1E-18,1E-18,1E-18])
        self.model.addCompartment("C", 1,True)
        # Creating the random initial concetration if it's out of condition
        if i_sp is None: 
            i_sp=np.random.rand(self.MpDf.shape[0])
        elif(len(i_sp) !=self.MpDf.shape[0]):
            i_sp=np.random.rand(self.MpDf.shape[0])
        
        # Creating the random initial rates if it's out of condition
        if rt is None:  
            rt=np.random.rand(self.MpDf.shape[1])
        elif (len(rt) !=self.MpDf.shape[1]):
           rt=np.random.rand(self.MpDf.shape[1])
           
        # Adding species to the model and treshold concetration
        for i in range(len(self.SpIdStrArray)):
            self.model.addSpecies(self.SpIdStrArray[i], compartment="C", initConcentration=i_sp[i],forceRegenerate=True)
        # Adding reactions and reactions rate constants
        for i in range(self.MpDf.shape[1]):
            self.model.addParameter("k"+str(i),rt[i],True)
            
            prod=[]
            for j in np.where(self.MpDf.iloc[:,i])[0]:
                if self.MpDf.iloc[j,i].astype(int) > 1:
                    prod.append(str(self.MpDf.iloc[j,i].astype(int))+self.MpDf.index[j])
                else:
                    prod.append(self.MpDf.index[j])
            reac=[]
            rate="k"+str(i)
            for j in np.where(self.MrDf.iloc[:,i])[0]:
                if self.MrDf.iloc[j,i].astype(int) > 1:
                    reac.append(str(self.MrDf.iloc[j,i].astype(int))+self.MrDf.index[j])
                else:
                    reac.append(self.MrDf.index[j])
                    
                if self.MrDf.iloc[j,i]>1:
                    for k in range(1,self.MrDf.iloc[j,i].astype(int)+1):
                        rate+=" * "+self.MrDf.index[j]
                else:
                    rate+=" * "+self.MrDf.index[j]
   
            self.model.addReaction("r"+str(i), reac, prod, rate)
        
        # Adding events so if conentrations of species are under cutoff, 
        # will have zero concentration 
        # for i in range(len(self.SpIdStrArray)):
        #     self.model.addEvent("e"+str(i),True,self.SpIdStrArray[i]+"<"+str(cutoff),True)
        #     self.model.addEventAssignment("e"+str(i),self.SpIdStrArray[i],"0",True)
        #     self.model.addTrigger("e"+str(i),self.SpIdStrArray[i]+"<"+str(cutoff),forceRegenerate=True)
        
        self.model.regenerateModel()
        self.IsSbmlbool=False
        
    # Function that load dynamical model directly form the sbml file     
    def setSmblModel(self):
        if not self.IsSbmlbool:
            print("network do not correspond to an sbml file")
            return
        # self.model = te.loadSBMLModel(self.fname)
        self.model = re.RoadRunner(self.fname)
    
        
    # Function that return concentration species in right order and cosidering a cutoff
    def getSpConArray(self,cutoff=None):
        
        if cutoff is None:
            cutoff=self.CutoffFloat
        con=[]
        for i in self.SpIdStrArray:
            con.append(self.model.getValue(i))
        con=np.array(con)
        con[con<cutoff]=0
        return con
        
    # Function that return reaction rates in right order and cosidering a cutoff
    def getRpRateArray(self,cutoff=0.1):
        con=[]
        for i in range(self.MrDf.shape[1]):
            con.append(self.model.getValue('r'+str(i)))
        con=np.array(con)
        con[con<cutoff]=0
        return 
        
    # returns the boolean abstract state (closure) for a given reaction network rn
    # species concentrations vector s and flow constants vector f
    def getSpAbstracArray(self,s,f,cutoff=None):
        if cutoff is None:
            cutoff=self.CutoffFloat
        k = self.getIndArrayFromBt(self.getClosureFromSp(sp_set=self.SpIdStrArray[np.where(s>=cutoff)[0]],r_set=np.where(f>0)[0],bt_type=True))
        v = np.repeat(False,len(s))
        v[k] = True
        return v
    
    # returns the boolean active reaction
    # species concentrations vector s
    def getActiveRpArray(self,s,cutoff=None):
        
        if cutoff is None:
            cutoff=self.CutoffFloat
        sp = np.where(s>cutoff)[0]        
        def ac_reac(reac):
            return all(np.isin(np.where(reac)[0],sp))
        
        k=np.where(np.apply_along_axis(ac_reac, 0, self.MrDf))[0]
        v = np.repeat(False,self.MrDf.shape[1])
        v[k] = True
        return v

    # returns the boolean active species 
    # species concentrations vector s
    def getActiveSpArray(self,s,cutoff=None):
        
        if cutoff is None:
            cutoff=self.CutoffFloat
        sp = np.where(s>cutoff)[0]
        def ac_reac(reac):
            return all(np.isin(np.where(reac)[0],sp))
        
        v = np.repeat(False,len(s))
        reac = np.where(np.apply_along_axis(ac_reac, 0, self.MrDf))[0]
        if len(reac)>0:
            sp = np.where((self.MrDf.iloc[:,reac].sum(axis=1)+self.MpDf.iloc[:,reac].sum(axis=1))>0)[0]
            v[sp] = True
            
        return v

    
    # Function that simulates the dynamics of the proposed model. It takes as 
    # input the initial time ti and the final time tf, and the number of steps 
    # to simulate by means of the variable steps, and a cutoff concentration  
    # that considers if a species is present or not. The function generates the 
    # dataframes RNIRG.con with RNIRG.rate, which correspond to concentrations and 
    # processes of the simulation respectively. It also generates tre more variables, 
    # RNIRG.abst, RNIRG.a_sp and RNIRG.a_r. Those corresponde to the abstraction, 
    # active species and active reaction. This dataframes will append another 
    # line if any of the species change their concentration below or over the cutoff  
    def runModel(self,ti=0,tf=50,steps=100,cutoff=None): 
        
        if cutoff is None:
            cutoff=self.CutoffFloat
        else:
            self.CutoffFloat=cutoff
        # Generation of kinetic contstant vector
        k=[]
        rate_name=[]
        for i in range(self.MrDf.shape[1]):
            k.append(self.model.getGlobalParameterByName('k'+str(i)))
            rate_name.append('r'+str(i))
        
        k=np.array(k)
        
        # Considering if the simulation has benn aready run, and generation of
        # output variables
        t_step=(tf-ti)/steps
        c_st=self.model.getFloatingSpeciesAmountsNamedArray()[0]
        sp_names=self.model.getFloatingSpeciesAmountsNamedArray().colnames
        c_st[c_st<cutoff]=0
        con=[[ti]+c_st.tolist()]
        reac=self.model.getReactionRates()
        ac_reac=self.getActiveRpArray(c_st,cutoff)
        reac[ac_reac==False]=0
        rate=[[ti]+reac.tolist()]
        
        # # Check if the abstraction has changed, if remain the same, no data is added.
        # try:
        #     if any(self.RN.abst[-1]!=(c_st>0)):
        #         abst=[[ti]+self.getSpAbstracArray(c_st,k,cutoff).tolist()]
        #         a_sp=[[ti]+self.getActiveSpArray(c_st,cutoff).tolist()]
        #         a_r=[[ti]+self.getActiveRpArray(c_st,cutoff).tolist()]
        # # Exception is generated if the RN.abst variable is not initialized beforehand.
        # except:
        #     abst=[[ti]+self.getSpAbstracArray(c_st,k,cutoff).tolist()]
        #     a_sp=[[ti]+self.getActiveSpArray(c_st,cutoff).tolist()]
        #     a_r=[[ti]+self.getActiveRpArray(c_st,cutoff).tolist()]
        
        # Running the simulation
        for i in range(steps):
            t=t_step*i+ti
            self.model.simulate(start=t,end=t+t_step,points=2)
            n_st=self.model.getFloatingSpeciesAmountsNamedArray()[0]
            n_st[n_st<cutoff]=0
            con.append([t+t_step]+n_st.tolist())
            reac=self.model.getReactionRates()
            
            ac_reac=self.getActiveRpArray(n_st,cutoff)
            reac[ac_reac==False]=0
            rate.append([t+t_step]+reac.tolist())
            # # appending data if the active species changes
            # if not all((c_st<cutoff)==(n_st<cutoff)):
            #     abst.append([t+t_step]+self.getSpAbstracArray(n_st,k,cutoff).tolist())
            #     a_sp.append([t+t_step]+self.getActiveSpArray(n_st,cutoff).tolist())
            #     a_r.append([t+t_step]+ac_reac.tolist())
            # c_st=n_st.copy()
            c_st[c_st<cutoff]=0
        
        
        # Generating output
        con=pd.DataFrame(con,columns=['time']+sp_names)
        con=con.set_index('time')
        
        # abst=pd.DataFrame(abst,columns=['time']+sp_names)
        # abst=abst.set_index('time')
        
        # a_sp=pd.DataFrame(a_sp,columns=['time']+sp_names)
        # a_sp=a_sp.set_index('time')
        
        
        rate=pd.DataFrame(rate,columns=['time']+rate_name)
        rate=rate.set_index('time')
        
        # a_r=pd.DataFrame(a_r,columns=['time']+rate_name)
        # a_r=a_r.set_index('time')
        
        
        try:
            self.SpConDf=pd.concat([self.SpConDf,con])
            self.RpRateDf=pd.concat([self.RpRateDf,rate])
            #self.abst=pd.concat([self.abst,abst])
            #self.a_sp=pd.concat([self.a_sp,a_sp])
            #self.a_r=pd.concat([self.a_r,a_r])
            
        except:
            self.SpConDf=con
            self.RpRateDf=rate
            # self.abst=abst
            # self.a_sp=a_sp
            # self.a_r=a_r
    
    def getAbstracDf(self,abst_type="non_null",cutoff=None):    
        
        # Setting cutoff
        if cutoff is None:
            cutoff=self.CutoffFloat
        elif cutoff < self.CutoffFloat:
            print("Input cutoff is lower that used on simulation, sumilation cutoff wil be used.")
            cutoff = self.CutoffFloat
        
        #  Initializaziation of variables
        if abst_type=="non_null":
            abst=pd.DataFrame()
            abst=pd.concat([abst,self.SpConDf.loc[self.SpConDf.index[0:1],:]>0])
        elif abst_type=="close_abst":
            abst=[]
            abst=[[self.SpConDf.index[0]]+self.getSpAbstracArray(self.SpConDf.iloc[0],self.RpRateDf.iloc[0],cutoff).tolist()]
        elif abst_type=="active_species":
            abst=[]
            abst=[[self.SpConDf.index[0]]+self.getActiveSpArray(self.SpConDf.iloc[0],cutoff).tolist()]
        elif abst_type=="active_reactions":
            abst=[]
            abst=[[self.SpConDf.index[0]]+self.getActiveRpArray(self.SpConDf.iloc[0],cutoff).tolist()]

        
        # for over the concentration dataframe
        for i in range(self.SpConDf.shape[0]-1):
            
      
                if abst_type=="non_null":
                    cs=self.SpConDf.iloc[i]
                    ns=self.SpConDf.iloc[i+1]
                    # condition of abstraction change the abstaraction are added
                    if any((cs>cutoff) != (ns>cutoff)) and i>0:
                        abst=pd.concat([abst,self.SpConDf.loc[self.SpConDf.index[i+1:i+2],:]>0])
                elif abst_type=="close_abst":
                    cs=self.getSpAbstracArray(self.SpConDf.iloc[i],cutoff)
                    ns=self.getSpAbstracArray(self.SpConDf.iloc[i+1],cutoff)
                    # condition of abstraction change the abstaraction are added
                    if any((cs>cutoff) != (ns>cutoff)) and i>0:
                        abst.append([self.SpConDf.index[i+1]]+self.getSpAbstracArray(self.SpConDf.iloc[i+1],self.RpRateDf.iloc[0],cutoff).tolist())
                elif abst_type=="active_species":
                    cs=self.getActiveSpArray(self.SpConDf.iloc[i],cutoff)
                    ns=self.getActiveSpArray(self.SpConDf.iloc[i+1],cutoff)
                    # condition of abstraction change the abstaraction are added
                    if any((cs>cutoff) != (ns>cutoff)) and i>0:
                        abst.append([self.SpConDf.index[i+1]]+self.getActiveSpArray(self.SpConDf.iloc[i+1],cutoff).tolist())
                elif abst_type=="active_reactions":
                    cs=self.getActiveRpArray(self.SpConDf.iloc[i],cutoff)
                    ns=self.getActiveRpArray(self.SpConDf.iloc[i+1],cutoff)
                    # condition of abstraction change the abstaraction are added
                    if any((cs>cutoff) != (ns>cutoff)) and i>0:
                        abst.append([self.SpConDf.index[i+1]]+self.getActiveRpArray(self.SpConDf.iloc[i+1],cutoff).tolist())
        
        rate_name=[]
        for i in range(self.MrDf.shape[1]):
            rate_name.append('r'+str(i))
        sp_names=self.model.getFloatingSpeciesAmountsNamedArray().colnames
        
        if abst_type!="non_null" and abst_type!="active_reactions":
            abst=pd.DataFrame(abst,columns=['time']+sp_names)
            abst=abst.set_index('time')
        elif abst_type!="non_null":
            abst=pd.DataFrame(abst,columns=['time']+rate_name)    
            abst=abst.set_index('time')
        return abst
    
    # Function that modify parameters of a given model. It recives as input a 
    # vector of the initial concentration of species i_sp, the reactive 
    # constant vector rt and the concentration 
    def setMakParam(self, i_sp=None ,rt=None,init=False):
        
        if self.IsSbmlbool:
            print("The initialized model corresponds to an sbml model, therefore the perturbation function cannot be used.")
            return
        
        # Creating the random initial concetration if it's out of condition
        if i_sp is None: 
            i_sp=np.random.rand(self.MpDf.shape[0])
        elif(len(i_sp) !=self.MpDf.shape[0]):
            i_sp=np.random.rand(self.MpDf.shape[0])
        
        # Creating the random initial concetration if it's out of condition
        if rt is None:  
            rt=np.random.rand(self.MpDf.shape[1])
        elif (len(rt) !=self.MpDf.shape[1]):
            rt=np.random.rand(self.MpDf.shape[1])
        
        if init:
            del(self.SpConDf)
            del(self.RpRateDf)
            # del(self.abst)
            # del(self.a_r)
            # del(self.a_sp)

            for i in range(self.MpDf.shape[1]):
                self.model.setGlobalParameterByName("k"+str(i),rt[i])
            
            for i in range(self.MpDf.shape[0]):
                self.model.setInitConcentration(self.model.getIds()[i], i_sp[i])
                
        else:
            
            for i in range(self.MpDf.shape[1]):
                self.model.setGlobalParameterByName("k"+str(i),rt[i])
            
            for i in range(self.MpDf.shape[0]):
                self.model.setValue(self.model.getIds()[i], i_sp[i])
    
    # randomize components of a state or flow vector using a log normal distribution
    # v vector, mask of affected components (defaults to active ones), mu and sigma parameters of the distribution
    def getRandomizePert(self,v,mask=None,mu=0,sigma=.5):
        v_out=v.copy()
        if mask is None:
            mask=(v_out!=0)
        mask = np.repeat(True,len(v_out)) & mask
        k = np.where(mask)[0]
        if (len(k)==0): 
            return(v_out)
        v_out[k] = np.exp(np.random.normal(mu,sigma,len(k)))
        return v_out
    
    # function that generates a random perturbation of eplison size with respect 
    # v vector, only considering the mask species.
    def getStatePert(self,v,mask=None,epsilon=1):
        if mask is None:
            mask=(v!=0)
        mask = np.repeat(True,len(v)) & mask
        v_p=np.random.normal(0,0.5,len(v))
        v_p[~mask]=0
        v_p = v + epsilon*(v_p / np.linalg.norm(v_p))
        v_p[v_p<0]=0
        return v_p
    
    # activation of components of v (changing values from 0 to 1 or from >0 to 0), mask (defaults to all components) 
    # n total active components required, up to p active components to be preserved (defaults to all active components)
    def getPertActivation(self,v,mask=None,n=None,p=None):
        v_out=v.copy()
        if mask is None:
            mask=np.repeat(True,len(v_out))
        if n is None:
            n = np.sum(mask)
        if p is None:
            p = np.sum(v_out[mask]>0)
        
        w = v_out[mask]
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
  
        v_out[mask] = w
        return v_out

    # perturb a vector by adding or substracting components and randomizing the result
    # v vector, d component delta, at least nmin active components, sigma random dispersion
    def getPertAddAndRandomize(self,v,d=1,nmin=5,sigma=.5):
        n = np.max([np.sum(v>0)+d,nmin])
        v_out = self.getPertActivation(v,n=n)
        v_out = self.getRandomizePert(v_out,sigma=sigma)
        
        return v_out

    # perturb a vector by adding or substracting components and randomizing the result
    # v vector respect at his current position, d component gamma, at least nmin active components, 
    # min_epsilon minimun random size perturmation, max_epsilon maximum random size pertubartion
    def getPertAddAndState(self,v,d=1,nmin=5,max_epsion=1,min_epsilon=0.1):
        n = np.max([np.sum(v>0)+d,nmin])
        v_out = self.getPertActivation(v,n=n)
        epsilon=np.random.uniform(min_epsilon,max_epsion)
        v_out = self.getStatePert(v_out,epsilon=epsilon)
        return v_out    

    # creates a structure to study evolutive paths of a reaction network rn
    # result: the reaction network rn, a random walk list rw
    # each rw[[i]] contains matrices f, s, p, c, a and r were each column stores
    # the starting state (f,s), the perturbed state p, its convergence (c), its abstraction (a) and used species (u)
    # for each perturbation and simulation step applied
    def resetRw(self):
      self.RwListDictDf=[]
      self.RwListDictDf.append(dict(sim=None,f=None,s=None,p=None,c=None,a=None,ac=None,u=None))
    
    
       
    
    # working script to generate a random reaction network and perturbation/simulation random walks
    # if e is provided new random walks or steps are added to e
    # if rn is provided that reaction network is used instead of generating a random one
    # the random walks to be created or completed are defined by range w, by default 1:10
    # by default the number of perturbation and simulation steps l is set to 10
    # cutoff is the concentration threshold for a species to be reactive and n is the number of simulation steps
    # trys are the number of perturbations done if the integrator have covergence problems.
    # sim_save if the simulation (con, rate, abst, a_sp and a_r) of each simulation is stored
    # fname is the json filename where the results are stored.
    # the result is available in global variable e (the value is also returned by the function)    
    # t the time elapsed in the dynamic simulation is stored in the random walk
    # f a dataframe of the flow vectors in the random walk
    # p a dataframe the perturbed states before the simulation starts in the random walk
    # s a dataframe of the states before the perturbation
    # c a dataframe of the convergent states in the random walk
    # a a dataframe of the abstractions of the convergent state (closure)
    # ac a dataframe of the active species of the convergent state (closure)
    # u a dataframe of the active species in the random walk (used species)
    # sim a list of dataframe whit the simulation data (con, rate, abst, a_sp and a_r dataframes)
    def setRw(self,w=range(10),l=10,cutoff=.1,rt=None,n=5000,trys=10,sim_save=True,fname="rand_walk.json"):
                
        try:
            self.RwListDictDf
        except:
            self.resetRw() # the starting structure to store the random walks to be generated
        
        if w==None or w==1 or w==0:
            w=range(0,1)
            
        for i in w:  # for each random walk
            
            if i>len(self.RwListDictDf)-1: # this is a new random walk
                    self.RwListDictDf.append(dict(f=None,s=None,p=None,c=None,a=None,u=None,t=None)) # matrices are created to store the steps in columns
            if self.RwListDictDf[i]['f'] is None: # this is a void random walk (0 steps)
                s = np.zeros(len(self.SpIdStrArray)) # the current state is zeroed
                if rt is None:
                    f = self.getRandomizePert(np.ones(self.MrDf.shape[1])) # the flow vector is randomized (each random walk has a different f)
                else:
                    f = rt
                
                if i==0:
                    try:
                        del(self.SpConDf)
                        del(self.RpRateDf)
                        # del(self.abst)
                        # del(self.a_sp)
                        # del(self.a_r)
                        self.setMakModel(i_sp=s ,rt=f ,cutoff=cutoff) # initialization of the model
                    except:
                        self.setMakModel(i_sp=s ,rt=f ,cutoff=cutoff) # initialization of the model
                else:
                    self.setMakParam(i_sp=s ,rt=f,init=True)
                    
            else: # this is a random walk with a number of steps already accumulated
              j = self.RwListDictDf[i]['c'].shape[1]  # the number of steps up to now
              s = self.RwListDictDf[i]['c'][:,j] # the exploration continues from the last convergence state in the random walk
              f = self.RwListDictDf[i]['f'][:,j] # the last flow vector is conserved
            
            fail=False
            for k in range(l):  # for each j step in random walk i
                print("walk: "+str(i+1)+", step: "+str(k+1))
                # the current state is stored in the random walk
                
                for m in range(trys):
                    # try:
                    if k==0:
                        s = self.getPertAddAndState(np.zeros(len(self.SpIdStrArray)),d=1,nmin=1)
                        c_st=pd.DataFrame(self.model.getFloatingSpeciesConcentrationsNamedArray(),
                                                     columns=self.model.getFloatingSpeciesConcentrationsNamedArray().colnames).T
                    else:
                        s = self.getPertAddAndState(np.array(self.SpConDf.iloc[-1]),d=1,nmin=1) # a delta perturbation is applied to current state
                        c_st = self.SpConDf.iloc[-1]
                    
                    # end perturbations, start simulation:
                    self.setMakParam(s)
                    start = time.time() 
                    self.runModel(n*k,n*(k+1),n,cutoff) # the perturbed state is simulated reaching a convergence state
                    end = time.time()
                    st=end-start 
                        # end simulation
                    break
                    
                    # except:
                    #     print("Convergence faliure, try number: ",m+1)
                
                iter_steps=k
                if (m+1)>=trys:
                    fail=True
                    break
                        
                if k==0:
                    self.RwListDictDf[i]['s']=c_st
                    self.RwListDictDf[i]['f'] = self.RpRateDf.iloc[-1].copy() # the flow vector is stored in the random walk
                    self.RwListDictDf[i]['p'] = pd.DataFrame(s,index=self.SpIdStrArray)  # the perturbed current state is stored in the random walk
                    self.RwListDictDf[i]['c'] = self.SpConDf.iloc[-1].copy() # the convergent state is stored in the random walk
                    self.RwListDictDf[i]['a'] = pd.DataFrame(self.getSpAbstracArray(self.SpConDf.iloc[-1],cutoff),
                                                             index=self.SpIdStrArray) # the abstraction is stored  
                    self.RwListDictDf[i]['ac'] = pd.DataFrame(self.getActiveSpArray(self.SpConDf.iloc[-1],cutoff),
                                                             index=self.SpIdStrArray) # a second abstraction is stored (active species)
                    self.RwListDictDf[i]['u'] = pd.DataFrame(s>0,index=self.SpIdStrArray)  # a thrid abstraction is stored (active species)
                    self.RwListDictDf[i]['t'] = [st] # the time elapsed in the dynamic simulation is stored in the random walk
                else:
                    self.RwListDictDf[i]['s'] = pd.concat([self.RwListDictDf[i]['s'],self.SpConDf.iloc[-1]],axis=1)
                    self.RwListDictDf[i]['f'] = pd.concat([self.RwListDictDf[i]['f'],self.RpRateDf.iloc[-1]],axis=1) # the flow vector is stored in the random walk
                    self.RwListDictDf[i]['p'] = pd.concat([self.RwListDictDf[i]['p'],pd.DataFrame(s,index=self.SpIdStrArray)],axis=1)  # the perturbed current state is stored in the random walk
                    self.RwListDictDf[i]['c'] = pd.concat([self.RwListDictDf[i]['c'],self.SpConDf.iloc[-1]],axis=1) # the convergent state is stored in the random walk
                    self.RwListDictDf[i]['a'] = pd.concat([self.RwListDictDf[i]['a'],pd.DataFrame(self.getSpAbstracArray(self.SpConDf.iloc[-1],cutoff),
                                                                                                  index=self.SpIdStrArray)],axis=1) # the abstraction is stored  
                    self.RwListDictDf[i]['ac'] = pd.concat([self.RwListDictDf[i]['ac'],pd.DataFrame(self.getActiveSpArray(self.SpConDf.iloc[-1],cutoff),
                                                                                                    index=self.SpIdStrArray)],axis=1) # a second abstraction is stored (active species) 
                    self.RwListDictDf[i]['u'] = pd.concat([self.RwListDictDf[i]['u'],pd.DataFrame(s>0,index=self.SpIdStrArray)],axis=1) # a thrid abstraction is stored (used species)
                    self.RwListDictDf[i]['t'].append(st) # the time elapsed in the dynamic simulation is stored in the random walk
                 
            
            if sim_save:
                h=0
                self.RwListDictDf[i]['sim']=[]
                for k in range(l):
                    self.RwListDictDf[i]['sim'].append(dict(con=self.SpConDf.iloc[n*(k+h):n*(k+h+1)],rate=self.RpRateDf.iloc[n*(k+h):n*(k+h+1)]))
                    h+=1
            
            if fail:
                if iter_steps != 0:
                    iter_steps-=1
                    
            self.RwListDictDf[i]['f'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['s'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['p'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['c'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['a'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['ac'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['u'].columns=range(iter_steps+1)
           
            
        out=copy.deepcopy(self.RwListDictDf)
        for i in out:
            
            i['s']=i['s'].to_json()
            i['f']=i['f'].to_json()
            i['p']=i['p'].to_json()
            i['c']=i['c'].to_json()
            i['a']=i['a'].to_json()
            i['ac']=i['ac'].to_json()
            i['u']=i['u'].to_json()
            for j in range(len(i['sim'])):
                for k in i['sim'][j]:
                    i['sim'][j][k].reset_index(inplace=True)
                    i['sim'][j][k]=i['sim'][j][k].to_json()
            
        with open(fname, "w") as outfile:
            json.dump(out, outfile)
        
        return
    
    # working script to generate a random reaction network and perturbation/simulation random walks
    # if e is provided new random walks or steps are added to e
    # if rn is provided that reaction network is used instead of generating a random one
    # the random walks to be created or completed are defined by range w, by default 1:10
    # by default the number of perturbation and simulation steps l is set to 10
    # cutoff is the concentration threshold for a species to be reactive and n is the number of simulation steps
    # trys are the number of perturbations done if the integrator have covergence problems.
    # sim_save if the simulation (con, rate, abst, a_sp and a_r) of each simulation is stored
    # fname is the json filename where the results are stored.
    # the result is available in global variable e (the value is also returned by the function)    
    # t the time elapsed in the dynamic simulation is stored in the random walk
    # f a dataframe of the flow vectors in the random walk
    # p a dataframe the perturbed states before the simulation starts in the random walk
    # s a dataframe of the states before the perturbation
    # c a dataframe of the convergent states in the random walk
    # a a dataframe of the abstractions of the convergent state (closure)
    # ac a dataframe of the active species of the convergent state (closure)
    # u a dataframe of the active species in the random walk (used species)
    # sim a list of dataframe whit the simulation data (con, rate, abst, a_sp and a_r dataframes)
    def setRwMulticore(self,w=range(10),l=10,cutoff=.1,rt=None,n=5000,trys=10,sim_save=True,fname="rand_walk.json"):
                
        try:
            self.RwListDictDf
        except:
            self.resetRw() # the starting structure to store the random walks to be generated
        
        if w==None or w==1 or w==0:
            w=range(0,1)
            
        for i in w:  # for each random walk
            
            if i>len(self.RwListDictDf)-1: # this is a new random walk
                    self.RwListDictDf.append(dict(f=None,s=None,p=None,c=None,a=None,u=None,t=None)) # matrices are created to store the steps in columns
            if self.RwListDictDf[i]['f'] is None: # this is a void random walk (0 steps)
                s = np.zeros(len(self.SpIdStrArray)) # the current state is zeroed
                if rt is None:
                    f = self.getRandomizePert(np.ones(self.MrDf.shape[1])) # the flow vector is randomized (each random walk has a different f)
                else:
                    f = rt
                
                if i==0:
                    try:
                        del(self.SpConDf)
                        del(self.RpRateDf)
                        # del(self.abst)
                        # del(self.a_sp)
                        # del(self.a_r)
                        self.setMakModel(i_sp=s ,rt=f ,cutoff=cutoff) # initialization of the model
                    except:
                        self.setMakModel(i_sp=s ,rt=f ,cutoff=cutoff) # initialization of the model
                else:
                    self.setMakParam(i_sp=s ,rt=f,init=True)
                    
            else: # this is a random walk with a number of steps already accumulated
              j = self.RwListDictDf[i]['c'].shape[1]  # the number of steps up to now
              s = self.RwListDictDf[i]['c'][:,j] # the exploration continues from the last convergence state in the random walk
              f = self.RwListDictDf[i]['f'][:,j] # the last flow vector is conserved
            
            fail=False
            for k in range(l):  # for each j step in random walk i
                print("walk: "+str(i+1)+", step: "+str(k+1))
                # the current state is stored in the random walk
                
                for m in range(trys):
                    # try:
                    if k==0:
                        s = self.getPertAddAndState(np.zeros(len(self.SpIdStrArray)),d=1,nmin=1)
                        c_st=pd.DataFrame(self.model.getFloatingSpeciesConcentrationsNamedArray(),
                                                     columns=self.model.getFloatingSpeciesConcentrationsNamedArray().colnames).T
                    else:
                        s = self.getPertAddAndState(np.array(self.SpConDf.iloc[-1]),d=1,nmin=1) # a delta perturbation is applied to current state
                        c_st = self.SpConDf.iloc[-1]
                    
                    # end perturbations, start simulation:
                    self.setMakParam(s)
                    start = time.time() 
                    self.runModel(n*k,n*(k+1),n,cutoff) # the perturbed state is simulated reaching a convergence state
                    end = time.time()
                    st=end-start 
                        # end simulation
                    break
                    
                    # except:
                    #     print("Convergence faliure, try number: ",m+1)
                
                iter_steps=k
                if (m+1)>=trys:
                    fail=True
                    break
                        
                if k==0:
                    self.RwListDictDf[i]['s']=c_st
                    self.RwListDictDf[i]['f'] = self.RpRateDf.iloc[-1].copy() # the flow vector is stored in the random walk
                    self.RwListDictDf[i]['p'] = pd.DataFrame(s,index=self.SpIdStrArray)  # the perturbed current state is stored in the random walk
                    self.RwListDictDf[i]['c'] = self.SpConDf.iloc[-1].copy() # the convergent state is stored in the random walk
                    self.RwListDictDf[i]['a'] = pd.DataFrame(self.getSpAbstracArray(self.SpConDf.iloc[-1],cutoff),
                                                             index=self.SpIdStrArray) # the abstraction is stored  
                    self.RwListDictDf[i]['ac'] = pd.DataFrame(self.getActiveSpArray(self.SpConDf.iloc[-1],cutoff),
                                                             index=self.SpIdStrArray) # a second abstraction is stored (active species)
                    self.RwListDictDf[i]['u'] = pd.DataFrame(s>0,index=self.SpIdStrArray)  # a thrid abstraction is stored (active species)
                    self.RwListDictDf[i]['t'] = [st] # the time elapsed in the dynamic simulation is stored in the random walk
                else:
                    self.RwListDictDf[i]['s'] = pd.concat([self.RwListDictDf[i]['s'],self.SpConDf.iloc[-1]],axis=1)
                    self.RwListDictDf[i]['f'] = pd.concat([self.RwListDictDf[i]['f'],self.RpRateDf.iloc[-1]],axis=1) # the flow vector is stored in the random walk
                    self.RwListDictDf[i]['p'] = pd.concat([self.RwListDictDf[i]['p'],pd.DataFrame(s,index=self.SpIdStrArray)],axis=1)  # the perturbed current state is stored in the random walk
                    self.RwListDictDf[i]['c'] = pd.concat([self.RwListDictDf[i]['c'],self.SpConDf.iloc[-1]],axis=1) # the convergent state is stored in the random walk
                    self.RwListDictDf[i]['a'] = pd.concat([self.RwListDictDf[i]['a'],pd.DataFrame(self.getSpAbstracArray(self.SpConDf.iloc[-1],cutoff),
                                                                                                  index=self.SpIdStrArray)],axis=1) # the abstraction is stored  
                    self.RwListDictDf[i]['ac'] = pd.concat([self.RwListDictDf[i]['ac'],pd.DataFrame(self.getActiveSpArray(self.SpConDf.iloc[-1],cutoff),
                                                                                                    index=self.SpIdStrArray)],axis=1) # a second abstraction is stored (active species) 
                    self.RwListDictDf[i]['u'] = pd.concat([self.RwListDictDf[i]['u'],pd.DataFrame(s>0,index=self.SpIdStrArray)],axis=1) # a thrid abstraction is stored (used species)
                    self.RwListDictDf[i]['t'].append(st) # the time elapsed in the dynamic simulation is stored in the random walk
                 
            
            if sim_save:
                h=0
                self.RwListDictDf[i]['sim']=[]
                for k in range(l):
                    self.RwListDictDf[i]['sim'].append(dict(con=self.SpConDf.iloc[n*(k+h):n*(k+h+1)],rate=self.RpRateDf.iloc[n*(k+h):n*(k+h+1)]))
                    h+=1
            
            if fail:
                if iter_steps != 0:
                    iter_steps-=1
                    
            self.RwListDictDf[i]['f'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['s'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['p'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['c'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['a'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['ac'].columns=range(iter_steps+1)
            self.RwListDictDf[i]['u'].columns=range(iter_steps+1)
           
            
        out=copy.deepcopy(self.RwListDictDf)
        for i in out:
            
            i['s']=i['s'].to_json()
            i['f']=i['f'].to_json()
            i['p']=i['p'].to_json()
            i['c']=i['c'].to_json()
            i['a']=i['a'].to_json()
            i['ac']=i['ac'].to_json()
            i['u']=i['u'].to_json()
            for j in range(len(i['sim'])):
                for k in i['sim'][j]:
                    i['sim'][j][k].reset_index(inplace=True)
                    i['sim'][j][k]=i['sim'][j][k].to_json()
            
        with open(fname, "w") as outfile:
            json.dump(out, outfile)
