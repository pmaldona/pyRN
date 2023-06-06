#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 12:43:52 2023

@author: abassi
"""

import numpy as np

# np.random.choice((1,2,3),p=(.2,0,.8))

import string
from scipy.stats import binom


class gen_rn_csp:

    def __init__( self,
                  bsp_N=10,bsp=None,bsp_wf=1/3,bsp_w=None,
                  sp_len_max=10, sp_paired_p=None,
                  sp2_wf=1/3, sp2_len_w=None, sp2_right_side_p=.5,
                  new_sp_len_w=None,
                  synt_reac_p=.5, ):
        if bsp is None: self.bsp = list(string.ascii_letters)[0:bsp_N]
        else: self.nsp = bsp
        if bsp_w is None:
            bsp_w = np.exp(np.arange(len(self.bsp))*np.log(bsp_wf))
        self.bsp_w = bsp_w/np.sum(bsp_w)
        if sp_paired_p is None:
            sp_paired_p = np.arange(sp_len_max,0,-1)/sp_len_max
        self.sp_paired_p = sp_paired_p
        l = sp_paired_p.shape[0]
        if sp2_len_w is None:
            sp2_len_w = np.exp(np.arange(l)*np.log(sp2_wf))
        self.sp2_len_w = sp2_len_w/np.sum(sp2_len_w)
        self.sp2_right_side_p = sp2_right_side_p
        if new_sp_len_w is None:
            new_sp_len_w = binom.pmf(list(range(l)),l-1,.5)
        self.new_sp_len_w = new_sp_len_w/sum(new_sp_len_w)
        self.synt_reac_p = synt_reac_p
        self.sp = []
        self.sp_len = []
        self.sp_nreac = []
        self.sp_nprod = []
        self.sp_count = 0
        self.spndx = {}
        self.reac = []
        
    def add_sp(self,l=None,s=None,reuse=False):
        if l is None and s is None:
            l = np.random.choice(self.new_sp_len_w.shape[0],p=self.new_sp_len_w) + 1
        i = None
        if s is None and reuse:
            k = np.argwhere(np.array(self.sp_len)==l)
            if k.shape[0]>0: i = np.random.choice(k[:,0])
        if i is None:
            if s is None:  
                s = ''.join(np.random.choice(self.bsp,size=l,replace=True,p=self.bsp_w))
            if s in self.spndx: i = self.spndx[s]
        if i is None:    
            i = self.sp_count
            self.spndx[s] = i
            self.sp.append(s)
            self.sp_len.append(len(s))
            self.sp_nreac.append(0)
            self.sp_nprod.append(0)
            self.sp_count = i + 1
        return i
        
    def rand_reaction(self,sp1=None,reac_p=None,reuse=False):
        if sp1 is None: sp1 = self.add_sp(reuse=reuse)            
        l1 = self.sp_len[sp1]
        l2 = np.random.choice(self.sp2_len_w.shape[0],p=self.sp2_len_w) + 1
        if l1>self.sp_paired_p.shape[0]: is_paired = False
        elif l1==1: is_paired = True
        else: is_paired = np.random.uniform()<self.sp_paired_p[l1-1]
        if reac_p is not None:
            if is_paired: synt_reac_p = reac_p
            else: synt_reac_p = 1 - reac_p
        else: synt_reac_p = self.synt_reac_p
        is_synt_reac = np.random.uniform()<synt_reac_p
        sp2_right_side = np.random.uniform()<self.sp2_right_side_p
        if is_paired:
            sp2 = self.add_sp(l2,reuse=reuse)
            if not sp2_right_side: sp1, sp2 = sp2, sp1
            sp3 = self.add_sp(s=self.sp[sp1]+self.sp[sp2])
        else:
            sp3 = sp1
            l2 = min(l2,l1-1)
            if sp2_right_side: l1 = l1 - l2
            else: l1 = l2
            sp1 = self.add_sp(s=self.sp[sp3][:l1])
            sp2 = self.add_sp(s=self.sp[sp3][l1:])
        if is_synt_reac: r = ((sp1,sp2),(sp3,))
        else: r = ((sp3,),(sp1,sp2))
        self.reac.append(r)
        for i in r[0]: self.sp_nreac[i] = self.sp_nreac[i] + 1
        for i in r[1]: self.sp_nprod[i] = self.sp_nprod[i] + 1
        return r
    
    def add_reactions(self,n=1,close_sp=False,reuse=False):
        if not close_sp:
            for i in range(n): self.rand_reaction(reuse=reuse)
        else:
            for i in range(n):
                k = np.argwhere(np.logical_or(np.array(self.sp_nreac)==0,np.array(self.sp_nprod)==0))
                if k.shape[0]==0: break
                j = np.random.choice(k[:,0])
                if self.sp_nreac[j]==0: self.rand_reaction(j,True,reuse)
                else: self.rand_reaction(j,False,reuse)

    def get_matrices(self):
        mr = np.zeros((len(self.sp),len(self.reac)))
        mp = np.zeros((len(self.sp),len(self.reac)))
        for k,r in enumerate(self.reac):
            for i in r[0]: mr[i,k] = mr[i,k] + 1
            for i in r[1]: mp[i,k] = mp[i,k] + 1
        return (mr,mp)
    
    def display(self):
        for r in self.reac:
            p = ""
            for i,s in enumerate(r[0]):
                if i>0: p = p + " + "
                p = p + self.sp[s]
            p = p + " => "
            for i,s in enumerate(r[1]):
                if i>0: p = p + " + "
                p = p + self.sp[s]
            print(p)

    def test(self):
        self.add_reactions(3)
        self.add_reactions(100,True)
        self.display()
        