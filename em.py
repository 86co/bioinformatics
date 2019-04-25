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

D = collections.namedtuple("D", "X, C")
GS_result = collections.namedtuple("GS_Result", "motif, score")    #Global search result

def main():
    Entrez.email = "fujioka86co@gmail.com"
    
    #sample set
    dna0=list("CGCAGTACGGCTATGTCAT")
    dna1=list("CGTAGCTGAGCTGCTAGCT")
    
    dataset=[]
    dataset.append(D(dna0,1))
    dataset.append(D(dna1,0))
    
    gmm = calc_gmm(dataset)
    theta_b = calc_fm(dataset)
    
    width = 5
    lmbd = 0.8
    
    b = 1.0
    
    motif, maxScore = global_search(dataset, gmm, theta_b, width, lmbd, b)
    print(motif, maxScore)
    
def calc_gmm(dataset):
    den=0
    num=0
    for d in dataset:
        den+=1
        if d.C==1: num+=1 
    gmm = num/den if den>0 else None
    return gmm

def calc_fm(dataset):
    num=[0,0,0,0]
    den=0
    for d in dataset:
        for x in d.X:
            num[to_int[x]]+=1
            den+=1
    theta_b = np.array([num[0], num[1], num[2], num[3]], float)/den if den>0 else None
    return theta_b

def global_search(dataset, gmm, theta_b, width, lmbd, b):
    gs_result=[]
    for d in dataset:
        for i in range(len(d.X)-width+1):
            motif = d.X[i:i+width]
            theta_m = to_psfm(motif, b)
            score = objective(lmbd, theta_m, gmm, theta_b, dataset)
            gs_result.append(GS_result(motif, score))
    motif_maxScore=""
    maxScore=0
    for res in gs_result:
        if res.score>maxScore:
            motif_maxScore=res.motif
            maxScore=res.score
    return motif_maxScore, maxScore
        
def to_psfm(motif, b):
    theta_m = np.zeros((len(SGM),len(motif)))
    for i, m_i in enumerate(motif):
        for a, sgm_a in enumerate(SGM):
            theta_m[a,i] = float(m_i==sgm_a)
    theta_m+=b/len(SGM)
    theta_m/=1+b
    return theta_m
        
def objective(lmbd, theta_m, gmm, theta_b, dataset):
    v = gmm*lmbd/(1-gmm*lmbd)
    q = gmm*(1-lmbd)/(1-gmm*lmbd)
    score = 0
    for d in dataset:
        mu = 0
        for i in range(len(d.X)-theta_m.shape[1]+1):
            prod=1
            for j in range(theta_m.shape[1]):
                prod*=theta_m[to_int[d.X[i+j]],j]/theta_b[to_int[d.X[i+j]]]
            mu+=prod
        mu/=len(d.X)-theta_m.shape[1]+1
        sig_y = 1/(1+1/v*mu)
        p = (1-q)*(1-sig_y) if d.C==0 else sig_y+q*(1-sig_y)
        score+=p
    return score
    
if __name__ == "__main__":
    main()