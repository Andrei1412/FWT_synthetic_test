#!/usr/bin/env python
# Generated with SMOP  0.41
#from libsmop import *
#import os
import numpy as np
#import matplotlib.pyplot as plt
#import pdb; 
def conv_wf(eii3_i,eiib_i,nt,dt):
    """
         Convolve forward and backward
    """
    (nt,ny,nz,nx) = eii3_i.shape
    sum_xx= np.zeros((ny,nz,nx))
    for it in range(1,nt):
#        sum_xx=sum_xx+dt*np.multiply(eii3[:,:,:,it-1],eiib[:,:,:,nt-it])
        sum_xx=sum_xx+dt*np.multiply(eii3_i[it-1,:,:,:],eiib_i[nt-it,:,:,:])     
    return sum_xx

def Vsp_read(nx,ny,nz,snap_file):
    fid=open(snap_file,'r')
    sek=np.fromfile(fid, dtype='<f4')
    print(sek.shape)
    Vp=np.reshape(sek,[ny,nz,nx])
    return Vp
    


# Read_Dev_strain_ex2.m
nx=200
ny=200
nz=75
nts=5000
nt=500
nxyz=nx*ny*nz

with open('../../Dev_Strain/fwd01_xyzts.eii3', 'rb') as fid:
    data_array = np.fromfile(fid, np.float32)
    data_array=data_array[0:14]
    dx=data_array[8]
    dy=data_array[9]
    dz=data_array[10]
    dt=data_array[11]
    del data_array
    
snap_file1='../../../Model/Models/vs3dfile_in.s'
snap_file2='../../../Model/Models/vp3dfile_in.p'
snap_file3='../../../Model/Models/rho3dfile_in.d'

Vs=Vsp_read(nx,ny,nz,snap_file1)
Vp=Vsp_read(nx,ny,nz,snap_file2)
rho=Vsp_read(nx,ny,nz,snap_file3)

Lamda=np.multiply(rho,(np.multiply(Vp,Vp)-2*np.multiply(Vs,Vs)));
Mu=np.multiply(rho,np.multiply(Vs,Vs))
Kappa=Lamda+2/3*Mu

#mainfolder_fig='/home/user/workspace/GMPlots/Figs/'
#eijp_names = {'eii3','exxp','eyyp','exyp','exzp','eyzp'};
eijp_names = {'exyp','exzp','eyzp'}
#eijp_names = {'eii3'}
fw_file='../../Dev_Strain/fwd01_xyzts.';
bw_file='../../Dev_Strain/adj01_xyzts.';

eij_name='eii3'
matrix_file_fw=fw_file+eij_name
matrix_file_bw=bw_file+eij_name
fid1=open(matrix_file_fw, 'rb')
matrix_dummy_fw=np.fromfile(fid1, np.float32)
matrix_dummy_fw=matrix_dummy_fw[15:]

fid2=open(matrix_file_bw, 'rb')
matrix_dummy_bw=np.fromfile(fid2, np.float32)  
matrix_dummy_bw=matrix_dummy_bw[15:]
eii3 = np.reshape(matrix_dummy_fw,[nt,ny,nz,nx]) 
eiib = np.reshape(matrix_dummy_bw,[nt,ny,nz,nx])      
  
sum1=9*conv_wf(eii3,eiib,nt,dt)
GK = - np.multiply(sum1,Kappa)   

for eij_name in eijp_names:
    matrix_file_fw=fw_file+eij_name
    matrix_file_bw=bw_file+eij_name
    fid1=open(matrix_file_fw, 'rb')
    matrix_dummy_fw=np.fromfile(fid1, np.float32)
    matrix_dummy_fw=matrix_dummy_fw[15:]
    
    fid2=open(matrix_file_bw, 'rb')
    matrix_dummy_bw=np.fromfile(fid2, np.float32)  
    matrix_dummy_bw=matrix_dummy_bw[15:]
    
#    if eij_name=='eii3':
#        eii3 = np.reshape(matrix_dummy_fw,[nt,ny,nz,nx]) 
#        eiib = np.reshape(matrix_dummy_bw,[nt,ny,nz,nx])        
#        sum1=9*conv_wf(eii3,eiib,nt,dt)
#        GK = - np.multiply(sum1,Kappa)   

#    if eij_name=='exxp':
#        exxp = np.reshape(matrix_dummy_fw,[nt,ny,nz,nx]) 
#        exxb = np.reshape(matrix_dummy_bw,[nt,ny,nz,nx])  
#        sumxx=conv_wf(exxp,exxb,nt,dt)
#        
#    if eij_name=='eyyp':
#        eyyp = np.reshape(matrix_dummy_fw,[nt,ny,nz,nx]) 
#        eyyb = np.reshape(matrix_dummy_bw,[nt,ny,nz,nx])  
#        sumyy=conv_wf(eyyp,eyyb,nt,dt)   
#        sumxxyy=conv_wf(exxp+eyyp,exxb+eyyb,nt,dt)
#        del exxp 
#        del exxb
#        del eyyp
#        del eyyb
        
    if eij_name=='exyp':
        exyp = np.reshape(matrix_dummy_fw,[nt,ny,nz,nx]) 
        exyb = np.reshape(matrix_dummy_bw,[nt,ny,nz,nx])  
        sumxy=4*conv_wf(exyp+eii3,exyb+eiib,nt,dt)
        del exyp 
        del exyb
        
    if eij_name=='exzp':
        exzp = np.reshape(matrix_dummy_fw,[nt,ny,nz,nx]) 
        exzb = np.reshape(matrix_dummy_bw,[nt,ny,nz,nx])  
        sumxz=4*conv_wf(exzp+eii3,exzb+eiib,nt,dt)    
        del exzp 
        del exzb        
        
    if eij_name=='eyzp':
        eyzp = np.reshape(matrix_dummy_fw,[nt,ny,nz,nx]) 
        eyzb = np.reshape(matrix_dummy_bw,[nt,ny,nz,nx])  
        sumyz=4*conv_wf(eyzp+eii3,eyzb+eiib,nt,dt)     
        del eyzp 
        del eyzb   
#        del eii3 
#        del eiib  
        
#sum2 = sumxx+sumyy+sumxxyy+2*(sumxy+sumxz+sumyz)        
sum2 = 0.5*(sumxy+sumxz+sumyz)        
GM = -2*np.multiply(sum2,Mu)

GS=2*(GM-4/3*np.multiply(np.divide(Mu,Kappa),GK))
GP=2*np.multiply(np.divide((Kappa+4/3*Mu),Kappa),GK)	        
     
GS_arr=GS.reshape(-1)
np.savetxt('KS.txt', GS_arr)

GP_arr=GP.reshape(-1)
np.savetxt('KP.txt', GP_arr)

#print('max GM_arr'+str((GM_arr[0])))
#print('max GK_arr'+str((GK_arr[0])))
