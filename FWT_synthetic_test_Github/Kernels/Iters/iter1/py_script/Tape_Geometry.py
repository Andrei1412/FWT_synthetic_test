#!/usr/bin/env python
# Generated with SMOP  0.41
#from libsmop import *
#import os
import numpy as np
#import math
#import scipy.ndimage as ndimage
#import matplotlib.pyplot as plt
#import pdb; 
def tape_xyz_matrix(nx,ny,nz,R_xyz):
    Tp=np.ones((ny,nz,nx))
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            for k in range(1,nz+1):
                
                if (i<R_xyz):
                    Tp[j-1,k-1,i-1]*=1-np.exp(-0.5*i/R_nf)
                    
                if ((nx-i-1)<R_xyz):
                    Tp[j-1,k-1,i-1]*=1-np.exp(-0.5*(nx-i-1)/R_nf)
                    
                if (j<R_xyz):
                    Tp[j-1,k-1,i-1]*=1-np.exp(-0.5*j/R_nf)
                    
                if ((ny-j-1)<R_xyz):
                    Tp[j-1,k-1,i-1]*=1-np.exp(-0.5*(ny-j-1)/R_nf)   
                                        
    return Tp
   
# Read_Dev_strain_ex2.m
nx=200
ny=200
nz=75

nxyz=nx*ny*nz

dx=1.0
dy=dx
dz=dx
dt=0.02

##########################################################

#R_nf=50
R_nf=25
Gra_M=np.zeros((ny,nz,nx))
Gra_K=np.zeros((ny,nz,nx))

print('Tape matrix from boundary at distance = '+str(R_nf))     
Tp=tape_xyz_matrix(nx,ny,nz,R_nf)
Tp_arr=Tp.reshape(-1)
tape_file='Dump/Tp_xyz.txt'
np.savetxt(tape_file, Tp_arr)

print('finish creating tape xyz')
