import os
import sys

from Bio import Entrez
from Bio import SeqIO
from Bio import pairwise2

import math
import random
import numpy as np

import collections
sys.setrecursionlimit(10000)    #再起深度の上限

A="A"
T="T"
G="G"
C="C"
SGM=[A,T,G,C]

to_int = {
    A: 0,
    T: 1,
    G: 2,
    C: 3
}

seed_prior_weight = 1.0
motif_width = 9

D = collections.namedtuple("D", "X, C")
GS_result = collections.namedtuple("GS_Result", "psfm, score")    #Global search result

class Params():
    def __init__(self, dataset, lmbd=0.8):
        self.theta_m = np.zeros((len(SGM),motif_width))
        self.theta_b = self.calc_fm(dataset)
        self.lmbd    = lmbd
        self.gmm     = self.calc_gmm(dataset)
        
    def calc_fm(self, dataset):
        num=[0,0,0,0]
        den=0
        for d in dataset:
            for x in d.X:
                num[to_int[x]]+=1
                den+=1
        theta_b = np.array([num[0], num[1], num[2], num[3]], float)/den if den>0 else None
        return theta_b
        
    def calc_gmm(self, dataset):
        den=0
        num=0
        for d in dataset:
            den+=1
            if d.C==1: num+=1 
        gmm = num/den if den>0 else None
        return gmm
            
def main():
#    Entrez.email = "fujioka86co@gmail.com"
    
    #sample set
    dna0=list("CGCAGTACGGCTATGTCAT")
    dna1=list("CGTAGCTGAGCTGCTAGCT")
    dna2=list("TACTGCGTACTGCGATCGC")
    dna3=list("GCGTATGCGTACTGATACG")
    
    dataset=[]
    dataset.append(D(dna0,1))
    dataset.append(D(dna1,1))
    dataset.append(D(dna2,0))
    dataset.append(D(dna3,0))
    
    params = Params(dataset)
    
    params.theta_m = global_search(dataset, params)
    print(params.theta_m)
    params.theta_m, maxScore = local_search(dataset, params)

def global_search(dataset, params):
    gs_result=[]
    for d in dataset:
        for i in range(len(d.X)-motif_width+1):
            motif = d.X[i:i+motif_width]
            params.theta_m = to_psfm(motif)
            score = objective(dataset, params)
            gs_result.append(GS_result(params.theta_m, score))
#            print(motif, score)
    theta_m_maxScore=""
    maxScore=0
    for res in gs_result:
        if res.score>maxScore:
            theta_m_maxScore=res.psfm
            maxScore=res.score
    return theta_m_maxScore
        
def to_psfm(motif):
    theta_m = np.zeros((len(SGM),motif_width))
    for i, m_i in enumerate(motif):
        for a, sgm_a in enumerate(SGM):
            theta_m[a,i] = float(m_i==sgm_a)
    theta_m += seed_prior_weight/len(SGM)
    theta_m /= 1+seed_prior_weight
    return theta_m
        
def objective(dataset, params):
    v = params.gmm*params.lmbd/(1-params.gmm*params.lmbd)
    q = params.gmm*(1-params.lmbd)/(1-params.gmm*params.lmbd)
    score = 0
    
    for d in dataset:
        mu = 0
        for i in range(len(d.X)-motif_width+1):
            prod=1
            for j in range(motif_width):
                prod*=params.theta_m[to_int[d.X[i+j]],j]/params.theta_b[to_int[d.X[i+j]]]
            mu+=prod
        mu/=len(d.X)-params.theta_m.shape[1]+1
        sig_y = 1/(1+1/v*mu)
        p = (1-q)*(1-sig_y) if d.C==0 else sig_y+q*(1-sig_y)
        score+=p
        
    return score

def local_search(dataset, params):
    w = np.zeros((len(SGM), motif_width))
    for i in range(motif_width):
        for a in range(len(SGM)):
            w[a,i] = math.log(params.theta_m[a,i]/params.theta_b[a])
    
#    print(w)
    
    return w, 0
    
if __name__ == "__main__":
    main()