#!/usr/bin/env python
# Generated with SMOP  0.41
#from libsmop import *
#import os
import numpy as np
import os
import time
#import math
#import matplotlib.pyplot as plt
#import pdb; 
def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

def read_source(source_file):
    
    with open(source_file, 'r') as f:
        lines = f.readlines()    
    line0=lines[0].split()
    nShot=int(line0[0])
    S=np.zeros((nShot,3))
    for i in range(1,nShot+1):
        line_i=lines[i].split()
        S[i-1,0]=line_i[0]
        S[i-1,1]=line_i[1]
        S[i-1,2]=line_i[2]
    
    return nShot, S

def write_par_obs_source_i(Si):
    fname='e3d_mysource_i_obs.par'
    #file_default='e3d_mysource.par'
    os.system('cp e3d_obs.par e3d_mysource_i_obs.par')
    fid=open(fname,'a')
    fid.write("%s\n" %('xsrc='+str(Si[0])))
    fid.write("%s\n" %('ysrc='+str(Si[1])))
    fid.write("%s\n" %('zsrc='+str(Si[2])))

# Read_Dev_strain_ex2.m
nx=200
ny=200
nz=75
nts=5000
nt=500

#dx=0.4
dx=1.0
dy=dx
dz=dx
#dt=0.02
dt=0.05
source_file='../../StatInfo/SOURCE.txt'
nShot, S = read_source(source_file)

os.system('rm OutBin/*.*')
os.system('rm -r ../../Kernels/Vel_ob/*')

mkdir_p('../../Kernels/Vel_ob/Vel_ob_i')
############################################################
    
for ishot in range(1,nShot+1):
    dir_ob='../../Kernels/Vel_ob/Vel_ob_'+str(ishot)
    mkdir_p(dir_ob)
        
    write_par_obs_source_i(S[ishot-1,:])
    #time.sleep(5)
    job_file1 = 'fdrun-mpi_mysource_i_obs.sl'
    os.system("sbatch %s" %job_file1)
    time.sleep(40)
    os.system('mv ../../Kernels/Vel_ob/Vel_ob_i/*.* ../../Kernels/Vel_ob/Vel_ob_'+str(ishot))
    #exit()    
print('Finish Generate Data')
#####################################
time.sleep(20)


    
