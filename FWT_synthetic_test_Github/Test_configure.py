#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 11:00:02 2019

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt
#import os

def write_model_Vsp_name(rho,Vs,Vp,fname):
    [ny,nz,nx]=Vs.shape
    #[nx,nz,ny]=Vs.shape
    
    model_file='rho3d'+fname+'.d'
    fid=open(model_file,'wb')
    sek1 = np.array(rho, dtype=np.float32)
    sek1.astype('float32').tofile(fid)
    
    model_file0='vs3d'+fname+'.s'
    model_file1='vp3d'+fname+'.p'

    fid00=open(model_file0,'wb')
    sek2 = np.array(Vs, dtype=np.float32)
    sek2.astype('float32').tofile(fid00)

    fid01=open(model_file1,'wb')
    sek3 = np.array(Vp, dtype=np.float32)
    sek3.astype('float32').tofile(fid01)
    
def Vsp_read(nx,ny,nz,snap_file):
    fid=open(snap_file,'r')
    sek=np.fromfile(fid, dtype='<f4')
    Vp=np.reshape(sek,[ny,nz,nx])
    return Vp    

def write_rec_source(S,flag):
    if flag==1:
        flag_type='source'
        filename1='SOURCE.txt' 
    if flag==0:
        flag_type='station'
        filename1='STATION.txt'
        statnames =  ['A1','A2','A3','A4','A5','A6', 'A7','B1','B2','B3','B4','B5','B6', 'B7','C1','C2','C3','C4','C5','C6', 'C7','D1','D2','D3','D4','D5','D6', 'D7','E1','E2','E3','E4','E5','E6', 'E7','F1','F2','F3','F4','F5','F6', 'F7','G1','G2','G3','G4','G5','G6', 'G7']
        #statnames =['A1' ,'A2' ,'A3' ,'A4','A5' ,'A6' ,'A7' ,'A8','A9' ,'A10' ,'A11' ,'A12','A13' ,'A14' ,'A15' ,'A16','A17' ,'A18' ,'A19' ,'B1' ,'B2' ,'B3' ,'B4','B5' ,'B6' ,'B7' ,'B8','B9' ,'B10' ,'B11' ,'B12','B13' ,'B14' ,'B15' ,'B16','B17' ,'B18' ,'B19','C1' ,'C2' ,'C3' ,'C4','C5' ,'C6' ,'C7' ,'C8','C9' ,'C10' ,'C11' ,'C12','C13' ,'C14' ,'C15' ,'C16','C17' ,'C18' ,'C19' ,'D1' ,'D2' ,'D3' ,'D4','D5' ,'D6' ,'D7' ,'D8','D9' ,'D10' ,'D11' ,'D12','D13' ,'D14' ,'D15' ,'D16','D17' ,'D18' ,'D19']
    
    Ns=len(S)
    print(filename1)
    fid = open(filename1,'w')
    fid.write("%4d\n" %(Ns))    
    count=0
    while count<Ns:
        if flag==1:
            ss=str(count)
            sc=flag_type+ss
        else:
            sc=statnames[count]
        fid.write("%4d%4d%4d %10s\n" %(S[count,0],S[count,1],S[count,2],sc))        
        count=count+1
            
########################
nx=200
ny=200
nz=75

#dx=0.4
dx=1.0
dy=dx
dz=dx

x=dx*np.arange(1,nx+1,1)-dx
y=dy*np.arange(1,ny+1,1)-dy
z=dz*np.arange(1,nz+1,1)-dz

fname_true='file'
fname_init='file_h'

Vp=5.0*np.ones((ny,nz,nx))
Vs=0.5*Vp

#Vp=4.0*np.ones((ny,nz,nx))
rho=2.7*np.ones((ny,nz,nx))
write_model_Vsp_name(rho,Vs,Vp,fname_init)
#exit()
snap_file1='vs3dfile_h.s'
snap_file2='vp3dfile_h.p'
snap_file3='rho3dfile_h.d'

Vs1=Vsp_read(nx,ny,nz,snap_file1)
Vp1=Vsp_read(nx,ny,nz,snap_file2)
rho1=Vsp_read(nx,ny,nz,snap_file3)

K=Vs1[25,:,:]
fig = plt.figure(figsize=(6, 4))
#plt.imshow(K)
#plt.show()
plt.plot(K)
#################################
ix0=100;iz0=25;iy0=100;
for ix in range(1,nx+1):
    for iz in range(1,nz+1):
        for iy in range(1,ny+1):
            if (np.sqrt((ix-ix0-1)**2+(iz-iz0-1)**2+(iy-iy0-1)**2))<15:
                Vp[iy,iz,ix]=2.5
         
Vs=0.5*Vp           

write_model_Vsp_name(rho,Vs,Vp,fname_true)

snap_file1='vs3dfile.s'
snap_file2='vp3dfile.p'
snap_file3='rho3dfile.d'

Vs1=Vsp_read(nx,ny,nz,snap_file1)
Vp1=Vsp_read(nx,ny,nz,snap_file2)
rho1=Vsp_read(nx,ny,nz,snap_file3)

K=Vs1[25,:,:]
fig = plt.figure(figsize=(6, 4))
#plt.imshow(K)
#plt.show()
plt.plot(K)

#os.system('mv *file* Model/Model')
#print('Move files to Model/Model')
##################################
sx=np.arange(50,175+1,50)
sy=np.arange(50,175+1,50)


nsx=len(sx)
nsy=len(sy)

Ns=nsx*nsy
S=np.zeros((Ns,3))

for i in np.arange(1,nsx+1):
    for j in np.arange(1,nsy+1):
            S[(j-1)*nsx+i-1,0]=sx[i-1]
            S[(j-1)*nsx+i-1,1]=sy[j-1]
            S[(j-1)*nsx+i-1,2]=25
        
#S=[Sy, Sx, Sz]
###################################
#rx=np.arange(25,190+1,25)
#ry=np.arange(25,190+1,25)
            
rx=np.arange(40,160+1,20)
ry=np.arange(40,160+1,20)            

nrx=len(rx)
nry=len(ry)

Nr=nrx*nry

R=np.zeros((Nr,3))
for i in np.arange(1,nrx+1):
    for j in np.arange(1,nry+1):
        R[(j-1)*nrx+i-1,0]=rx[i-1]
        R[(j-1)*nrx+i-1,1]=ry[j-1]
        R[(j-1)*nrx+i-1,2]=1
        
        
write_rec_source(R,0)
write_rec_source(S,1)
##########################################
#Draw test configuration
plt.figure(figsize=(10,10))
for i in range(1,Ns+1):
    plt.plot((S[i-1,0])*dx,(S[i-1,1])*dx,'xr')
    
for j in range(1,Nr+1):
    plt.plot((R[j-1,0])*dx,(R[j-1,1])*dx,'ob')  
plt.title('Test configuration', loc='center') 
plt.xlim([0,nx*dx])
plt.ylim([0,ny*dy])
plt.grid()
plt.xlabel('x-axis (km)')
plt.ylabel('y-axis (km)')
#plt.gca().legend(('Observed','Simulated'))
plt.text(0, -15, 'x : sources at z='+str(S[0,2]*dx)+'km')
plt.text(40, -15, 'o : stations at z='+str(R[0,2]*dx-dx)+'km')
