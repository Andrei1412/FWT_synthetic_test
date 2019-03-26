#!/usr/bin/env python
# Generated with SMOP  0.41
#from libsmop import *
#import os
import numpy as np
#import math
import scipy.ndimage as ndimage
#import matplotlib.pyplot as plt
#import pdb; 
def read_source(source_file):
    
    with open(source_file, 'r') as f:
        lines = f.readlines()    
    line0=lines[0].split()
    nShot=int(line0[0])
    S=np.zeros((nShot,3))
    for i in range(1,nShot+1):
        line_i=lines[i].split()
        S[i-1,0]=int(line_i[0])
        S[i-1,1]=int(line_i[1])
        S[i-1,2]=int(line_i[2])
    
    return nShot, S


nx=200
ny=200
nz=75
nts=5000
nt=500
nxyz=nx*ny*nz

#dx=0.4
dx=1.0
dy=dx
dz=dx
dt=0.02

source_file='../../../StatInfo/SOURCE.txt'
station_file='../../../StatInfo/STATION.txt'

nShot, S = read_source(source_file)
nRec, R = read_source(station_file)                       
#with open(station_file, 'r') as f:
#    lines = f.readlines()
#    R=np.zeros((25,3))
#    i=0
#    for line in lines[1:]:
#        R1=line.split()
#        #S1.append([float(val) for val in line.split()])
#        for v in range(1,4):
#             R[i,v-1]=int(R1[v-1])
#        i=i+1
##########################################################

Gra_S=np.zeros((ny,nz,nx))
Gra_P=np.zeros((ny,nz,nx))

for ii in range(1,nShot+1):
    print('Sum Gradient Shot '+str(ii))  
    GS_file='All_shots/GS_shot'+str(ii)+'.txt'    
    GP_file='All_shots/GP_shot'+str(ii)+'.txt'
    
    GS_arr=np.loadtxt(GS_file)
    GS=np.reshape(GS_arr,[ny,nz,nx])
    GP_arr=np.loadtxt(GP_file)
    GP=np.reshape(GP_arr,[ny,nz,nx])
    
    #Tp=precondition_matrix(nx,ny,nz,S[ii-1,:],R,R_nf)
    tape_file='Dump/Tp_shot_'+str(ii)+'.txt'
    Tp_arr=np.loadtxt(tape_file)
    Tp=np.reshape(Tp_arr,[ny,nz,nx])
    
    Gra_S=Gra_S+np.multiply(GS,Tp)
    Gra_P=Gra_P+np.multiply(GP,Tp)

tape_file='Dump/Tp_xyz.txt'
Tp_arr=np.loadtxt(tape_file)
Tp=np.reshape(Tp_arr,[ny,nz,nx]) 

Gra_S=np.multiply(Gra_S,Tp)
Gra_P=np.multiply(Gra_P,Tp)

Gra_S[:,0:4,:]=0
Gra_P[:,0:4,:]=0

Gra_S=ndimage.gaussian_filter(Gra_S, 3)
Gra_P=ndimage.gaussian_filter(Gra_P, 3)

Gra_S_arr=Gra_S.reshape(-1)
np.savetxt('Gra_S.txt', Gra_S_arr)
    
Gra_P_arr=Gra_P.reshape(-1)
np.savetxt('Gra_P.txt', Gra_P_arr)
print('finish Summing Gradients')

print('max Gs='+str(np.max(Gra_S_arr)))
print('max Gp='+str(np.max(Gra_P_arr)))

