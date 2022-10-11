#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 10:14:52 2022

@author: pmaldona
"""

from .RNDS import RNDS
from .RNSRW import RNSRW
import numpy as np
from bitarray import bitarray as bt

class RNLI(RNSRW,RNDS):
    
    
    def setLogisticRN(self,reduc_rate=0.1):
        #set of overproduced species
        all_reac=bt(self.MrDf.shape[1])
        all_reac.setall(1)
        op_ind=self.getIndArrayFromBt(self.getOpSpBt(self.SpIdStrArray,all_reac))
        mr=self.MrDf[range(self.MrDf.shape[1])]
        
        # obtaining the max kinteic constant for each overproducible species
        for i in op_ind:
            op_sp_stoi=mr.iloc[i].max()
            op_max_reacs=np.where(mr.iloc[i]==op_sp_stoi)[0]
            max_kin=0
            if len(op_max_reacs)==1:
                max_kin=self.model.getValue("k"+str(op_max_reacs[0]))
            else:
                for j in op_max_reacs:
                    if self.model.getValue("k"+str(j))>max_kin:
                        max_kin=self.model.getValue("k"+str(j))

            # support and productive set update
            reac=bt(self.MrDf.shape[0])
            reac.setall(0)
            self.ProdListBt.append(reac)
            reac[i]=1
            self.ReacListBt.append(reac)
            
            # Stoichiometric reactive and productive part update.
            reac_v=np.zeros(self.MrDf.shape[0])
            self.MpDf[self.MpDf.shape[1]]=reac_v
            reac_v[i]=op_sp_stoi+1
            self.MrDf[self.MrDf.shape[1]]=reac_v
            
            #adding new parameters to the model
            self.model.addParameter("k"+str(self.MrDf.shape[1]-1),reduc_rate*max_kin,True)
            rate="k"+str(self.MrDf.shape[1]-1)
            
            for k in range(1,int(op_sp_stoi+2)):
                    rate+=" * "+self.SpIdStrArray[i]
            self.model.addReaction("r"+str(self.MrDf.shape[1]-1), [self.SpIdStrArray[i]], [], rate, True)
            
    @classmethod        
    def setFromTextWithMak(cls,file,SpConFileNameStr=None ,KConstFileNameStr=None ,cutoff=0.1,reduc_rate=0.1,log_fact=True):
        out=cls.setFromText(file)
        out.setMakModelFromFile(SpConFileNameStr ,KConstFileNameStr ,cutoff)
        if log_fact:
            out.setLogisticRN(reduc_rate)
        return(out)
    
    @classmethod
    def setFromSbmlWithMak(cls,file,modifiers=True,bond_con=True,rand_bc=False,SpConFileNameStr=None ,KConstFileNameStr=None ,cutoff=0.1,reduc_rate=0.1,log_fact=True):
        out=cls.setFromSbml(file,modifiers,bond_con,rand_bc)
        out.setMakModelFromFile(SpConFileNameStr ,KConstFileNameStr ,cutoff)
        if log_fact:
            out.setLogisticRN(reduc_rate)
        return(out)
    
    @classmethod
    def setRandomgeneratedWithInflowWithMak(cls,Nr=12,Ns=None,extra=.4,dist=lambda x: x*0+1,pr=0,pp=None,inflow=0.1,outflow=0.1,cutoff=0.1,
                                            filename=None, reduc_rate=0.1, log_fact=True):
        out=cls.setRandomgeneratedWithInflow(Nr,Ns,extra, dist, pr, pp, inflow, outflow)
        if filename is None:
            filename="Rand_net.txt"
        out.saveToText(filename)
        out.setMakModelFromFile(None ,None ,cutoff)
        if log_fact:
            out.setLogisticRN(reduc_rate)
        return(out)
   
    @classmethod
    def setRandomgeneratedNoInflowWithMak(cls,Nr=12,Ns=None,extra=.4,dist=lambda x: x*0+1,pr=0,pp=None,cutoff=0.1,
                                          filename=None,reduc_rate=0.1,log_fact=True):
        out=cls.setRandomgeneratedNoInflow(Nr,Ns,extra,dist,pr,pp)
        if filename is None:
            filename="Rand_net.txt"
        out.saveToText(filename)
        out.setMakModelFromFile(None ,None ,cutoff)
        if log_fact:
            out.setLogisticRN(reduc_rate)
        return(out)
    