#!/usr/bin/env python
# Generated with SMOP  0.41
#from libsmop import *
#import os
import numpy as np
import os
import time
from numpy.linalg import solve
#import math
#import matplotlib.pyplot as plt
#import pdb; 
def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

def write_par_source_i(Si):
    sname='FWD/e3d_mysource_i.par'
    #file_default='e3d_mysource.par'
    os.system('cp FWD/e3d_mysource_xyz_default.par FWD/e3d_mysource_i.par')
    fid=open(sname,'a')
    fid.write("%s\n" %('xsrc='+str(Si[0])))
    fid.write("%s\n" %('ysrc='+str(Si[1])))
    fid.write("%s\n" %('zsrc='+str(Si[2])))

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

def Vsp_read(nx,ny,nz,snap_file):
    fid=open(snap_file,'r')
    sek=np.fromfile(fid, dtype='<f4')
    Vp=np.reshape(sek,[ny,nz,nx])
    return Vp

def write_model_Vsp(rho,Vs,Vp):
    [ny,nz,nx]=Vs.shape
    
    model_file='rho3dfile1.d'
    fid=open(model_file,'wb')
    sek1 = np.array(rho, dtype=np.float32)
    sek1.astype('float32').tofile(fid)
    
    model_file0='vs3dfile1.s'
    model_file1='vp3dfile1.p'

    fid00=open(model_file0,'wb')
    sek2 = np.array(Vs, dtype=np.float32)
    sek2.astype('float32').tofile(fid00)

    fid01=open(model_file1,'wb')
    sek3 = np.array(Vp, dtype=np.float32)
    sek3.astype('float32').tofile(fid01)

def write_par_opt_source_i(Si):
    fname='FWD/e3d_mysource_i_opt.par'
    #file_default='e3d_mysource.par'
    os.system('cp FWD/e3d_opt_default.par FWD/e3d_mysource_i_opt.par')
    fid=open(fname,'a')
    fid.write("%s\n" %('xsrc='+str(Si[0])))
    fid.write("%s\n" %('ysrc='+str(Si[1])))
    fid.write("%s\n" %('zsrc='+str(Si[2])))

def Step_length_3D_pad2(Vs,Vp,rho,dx,dy,dz,dt,nt,Gra_Vp, Gra_Vs,Err,S,nShot):
    
    L20=Err
    L21=1.1*Err
    L22=0	
    flag=0
    Mu20=0
    Mu21=0
    Mu22=0	
    
    nx,nz,ny = Vs.shape
    #Step length
    st=0.1
    count=0
    #fact=2
    iter_max=8
    
    while (L21>L20 and count<iter_max ):
   
        if (count>0):
            #st=st/fact
            st=st/2	
        st_Vp=st*np.max(np.abs(Vp))/np.max(np.abs(Gra_Vp))
        st_Vs=st*np.max(np.abs(Vs))/np.max(np.abs(Gra_Vs))
                
        Vp0=Vp-st_Vp*Gra_Vp
        Vs0=Vs-st_Vs*Gra_Vs
        print('Mu21='+str(st))	
        print('max st_Vp0*Gra_Vp ='+str(st_Vp*np.max(np.abs(Gra_Vp))))  
        write_model_Vsp(rho,Vs0,Vp0)
        os.system('cp vs3dfile1.s Model/vs3dfile_opt.s')
        os.system('cp vp3dfile1.p Model/vp3dfile_opt.p')
        os.system('cp rho3dfile1.d Model/rho3dfile_opt.d')
        time.sleep(10)
        #
        #job_file = 'fdrun-mpi_mysource_opt.sl'
        for ishot in range(1,nShot+1): 
            write_par_opt_source_i(S[ishot-1,:])
            job_file = 'fdrun-mpi_ishot_opt.sl'
            os.system("sbatch %s" %job_file)
            time.sleep(40)
            os.system('mv ../../Vel_opt/Vel_ob_i/*.* ../../Vel_opt/Vel_ob_'+str(ishot))

        #job_file = '../../Read_GP_opt_calc_nShot.py'        
        job_file = '../../Read_GP_opt_calc_err_rwm_new.py'
        os.system("python %s" %job_file)
        time.sleep(70)
        f_err=open('err_opt.dat','r')
        Err1=np.fromfile(f_err,dtype='<f4')
        L21=Err1[1]
        Mu21=st
        print('L21')	
        print(L21)
        count=count+1
        if (count==iter_max):
            flag=1
            opt_step_L2=0
    ###
    count=0
    while (L21>L22 and flag==0):
        #st1=st*(1+1/fact)	 
        st=st*1.5
        st_Vp1=st*np.max(np.abs(Vp))/np.max(np.abs(Gra_Vp))
        st_Vs1=st*np.max(np.abs(Vs))/np.max(np.abs(Gra_Vs))
                
        Vp1=Vp-st_Vp1*Gra_Vp
        Vs1=Vs-st_Vs1*Gra_Vs
        print('Mu22='+str(st))	
        print('max st_Vp1*Gra_Vp='+str(st_Vp1*np.max(np.abs(Gra_Vp))))  
        write_model_Vsp(rho,Vs1,Vp1)
        os.system('cp vs3dfile1.s Model/vs3dfile_opt.s')
        os.system('cp vp3dfile1.p Model/vp3dfile_opt.p')
        os.system('cp rho3dfile1.d Model/rho3dfile_opt.d')
        time.sleep(10)        
        #
        for ishot in range(1,nShot+1): 
            write_par_opt_source_i(S[ishot-1,:])
            job_file = 'fdrun-mpi_ishot_opt.sl'
            os.system("sbatch %s" %job_file)
            time.sleep(40)
            os.system('mv ../../Vel_opt/Vel_ob_i/*.* ../../Vel_opt/Vel_ob_'+str(ishot))
            
        #job_file = '../../Read_GP_opt_calc_nShot.py'        
        job_file = '../../Read_GP_opt_calc_err_rwm_new.py'
        os.system("python %s" %job_file)
        time.sleep(70)
        f_err=open('err_opt.dat','r')
        Err1=np.fromfile(f_err,dtype='<f4')
        L22=Err1[1]
        print('L22')
        print(L22)
        Mu22=st
        count=count+1
        if (count==iter_max):
            flag=2
            opt_step_L2=st
        
    if (flag==0):
        A = [[Mu20**2, Mu20, 1], [Mu21**2, Mu21, 1], [Mu22**2, Mu22, 1]]    
        B = [L20,L21,L22]
        print(A)
        print(B)
        X=solve(A,B)

        opt_step_L2 = -X[1]/(2*X[0])
        if (opt_step_L2<0):
            opt_step_L2=Mu20
        if (opt_step_L2>Mu22):
            opt_step_L2=Mu22
    print('flag=')
    print(flag)
    return opt_step_L2,flag    
# Read_Dev_strain_ex2.m
nx=200
ny=200
nz=75
nts=5000
nt=500

dx=0.4
dy=0.4
dz=0.4
dt=0.02

iters=5
#nShot=16
source_file='../../../StatInfo/SOURCE.txt'
nShot, S = read_source(source_file)

os.system('rm -r ../../Vel_es/*')
os.system('rm -r ../../Vel_adj/*')
os.system('rm -r ../../Vel_opt/*')
os.system('rm -r /Dev_Strain/FW_xyz/*')

mkdir_p('../../Vel_es/Vel_es_i')
mkdir_p('../../Vel_adj/Vel_adj_i')
mkdir_p('../../Vel_opt/Vel_ob_i')
mkdir_p('../../Dev_Strain/FW_xyz')
mkdir_p('Dump')

############################################################
for it in range(2,iters+1):

    print('iNumber='+str(it))
    if (it==1):
        os.system('cp ../../../Model/Models/vs3dfile_h.s ../../../Model/Models/vs3dfile_in.s')
        os.system('cp ../../../Model/Models/vp3dfile_h.p ../../../Model/Models/vp3dfile_in.p')
        os.system('cp ../../../Model/Models/rho3dfile_h.d ../../../Model/Models/rho3dfile_in.d')
        
    if (it>1):
        os.system('cp %s' %'Dump/vs3dfile_iter_'+str(it-1)+'.s ../../../Model/Models/vs3dfile_in.s')
        os.system('cp %s' %'Dump/vp3dfile_iter_'+str(it-1)+'.p ../../../Model/Models/vp3dfile_in.p')
        os.system('cp %s' %'Dump/rho3dfile_iter_'+str(it-1)+'.d ../../../Model/Models/rho3dfile_in.d')
        
    time.sleep(5) 

    for ishot in range(1,nShot+1):
        print('isource='+str(ishot))
        if (it>0):
            dir_es='../../Vel_es/Vel_es_'+str(ishot)
            mkdir_p(dir_es)
            dir_adj='../../Vel_adj/Vel_adj_'+str(ishot)
            mkdir_p(dir_adj)            
            dir_opt='../../Vel_opt/Vel_ob_'+str(ishot)
            mkdir_p(dir_opt)
            dir_fwxyz='../../Dev_Strain/FW_xyz/FW_xyz_'+str(ishot)
            mkdir_p(dir_fwxyz)

        write_par_source_i(S[ishot-1,:])
        
        job_file11 = 'FWT_emod3d_shot_i_part1.sl'
        os.system("sbatch %s" %job_file11)
        time.sleep(130)
        os.system('mv ../../Vel_es/Vel_es_i/*.* ../../Vel_es/Vel_es_'+str(ishot))          
        os.system('mv ../../Dev_Strain/fwd01_xyzts.* ../../Dev_Strain/FW_xyz/FW_xyz_'+str(ishot)) 
        
    input('Copy adjoint sources and Press 1 to continue')
    
    for ishot in range(1,nShot+1):
        os.system('rm ../../../AdjSims/V3.0.7-a2a_xyz/Adj-InputAscii/*.*')
        os.system('cp ../../Vel_adj/Vel_adj_'+str(ishot)+'/*.* ../../../AdjSims/V3.0.7-a2a_xyz/Adj-InputAscii/')         
        fi=open('iShot.dat','w')
        (np.int64(ishot)).tofile(fi)
        fi.close() 
        print('isource='+str(ishot))
        
        os.system('cp ../../Dev_Strain/FW_xyz/FW_xyz_'+str(ishot)+'/fwd01_xyzts.* ../../Dev_Strain/')
        job_file12 = 'FWT_emod3d_shot_i_part2.sl'
        os.system("sbatch %s" %job_file12)
        time.sleep(130)

        job_file2 = 'kernel_shot_i_iter1.sl'
        os.system("sbatch %s" %job_file2)
        time.sleep(120)
        #os.system('mv ../../Vel_es/Vel_es_i/*.* ../../Vel_es/Vel_es_'+str(ishot))          
        os.system('mv KS.txt All_shots/GS_shot'+str(ishot)+'.txt')
        os.system('mv KP.txt All_shots/GP_shot'+str(ishot)+'.txt')
	#exit()
    fi1=open('iTape.dat','r')
    R_nf=np.int64(np.fromfile(fi1,dtype='int64'))
    fi1.close()
    print('R_nf='+str(R_nf))    
    #input('Press 2 and Enter')
    #os.system('python Sum_Gradients_nShot_vsp.py')

    job_file3 = 'sum_kernel.sl'
    os.system("sbatch %s" %job_file3)
    time.sleep(250)
    print('finish summing kernels')
    #os.system('python Model_update_CG_opt_step_iters_new.py')
    #################################
    snap_file1='../../../Model/Models/vs3dfile_in.s'
    snap_file2='../../../Model/Models/vp3dfile_in.p'
    snap_file3='../../../Model/Models/rho3dfile_in.d'
    Vs=Vsp_read(nx,ny,nz,snap_file1)
    Vp=Vsp_read(nx,ny,nz,snap_file2)
    rho=Vsp_read(nx,ny,nz,snap_file3)
       
    Gra_S_arr=np.loadtxt('Gra_S.txt')
#    Gra_S_arr=Gra_S_arr-np.mean(Gra_S_arr)
    Gra_Vs=np.reshape(Gra_S_arr,[ny,nz,nx])
    
    Gra_P_arr=np.loadtxt('Gra_P.txt')
#    Gra_P_arr=Gra_P_arr-np.mean(Gra_P_arr)    
    Gra_Vp=np.reshape(Gra_P_arr,[ny,nz,nx])

    
#    gra=np.zeros(2*nx*ny*nz)
#    grad_last=np.zeros(2*nx*ny*nz)
#    
#    gra[0:nx*ny*nz] = Gra_S_arr
#    gra[nx*ny*nz:2*nx*ny*nz] = Gra_P_arr
#    
#    if (it>0):
#        gra1=gra
#    
##    if (it>1):#CG Polak-Ribiere
##        grads_last_arr=np.loadtxt('Dump/Grads_iter_'+str(it-1)+'.txt')
##        gradp_last_arr=np.loadtxt('Dump/Gradp_iter_'+str(it-1)+'.txt')  
##        
##        grad_last[0:nx*ny*nz] = grads_last_arr
##        grad_last[nx*ny*nz:2*nx*ny*nz] = gradp_last_arr
##
##        beta_N=np.dot(gra,gra-grad_last)
##        beta_D=np.dot(gra,grad_last)
##        print('beta_N='+str(beta_N)+'and beta_D='+str(beta_D))
##        gra1=gra+beta_N/beta_D*grad_last        

    
    print('write gradient iter'+str(it))
    np.savetxt('Dump/Grads_iter_'+str(it)+'.txt', Gra_S_arr)   
    np.savetxt('Dump/Gradp_iter_'+str(it)+'.txt', Gra_P_arr)    
    ##############################
    #job_file = '../../Read_GP_obs_calc_nShot.py'        
    job_file = '../../Read_GP_obs_calc_err_rwm_new.py'
    os.system("python %s" %job_file)
    time.sleep(70) 
    f_err0=open('err_obs.dat','r')
    Err=np.fromfile(f_err0,dtype='<f4')
    print('Err='+str(Err[1]))
         
    os.system('cp err_obs.dat Dump/err_iter_'+str(it)+'.dat')   
    ##############################    
    opt_step_L2,flag = Step_length_3D_pad2(Vs,Vp,rho,dx,dy,dz,dt,nt,Gra_Vp, Gra_Vs,Err[1],S,nShot)
    print('opt_step_L2='+str(opt_step_L2))
    print('flag='+str(flag))
    if flag==1:
        print('No step length found!')
        break

    st_Vp3=opt_step_L2*np.max(np.abs(Vp))/np.max(np.abs(Gra_Vp))
    st_Vs3=opt_step_L2*np.max(np.abs(Vs))/np.max(np.abs(Gra_Vs))
            
    Vp3=Vp-st_Vp3*Gra_Vp
    Vs3=Vs-st_Vs3*Gra_Vs
    
    Vs_arr=Vs3.reshape(-1)
    Vs_max=6.0
    Vs_min=0.3
    for i in range(1,len(Vs_arr)+1):
        if Vs_arr[i-1]>Vs_max:
            Vs_arr[i-1]=Vs_max
        if Vs_arr[i-1]<Vs_min:
            Vs_arr[i-1]=Vs_min
    Vs2=np.reshape(Vs_arr,[ny,nz,nx])
    
    Vp_arr=Vp3.reshape(-1)
    Vp_max=10.0
    Vp_min=0.0
    for i in range(1,len(Vp_arr)+1):
        if Vp_arr[i-1]>Vp_max:
            Vp_arr[i-1]=Vp_max
        if Vp_arr[i-1]<Vp_min:
            Vp_arr[i-1]=Vp_min
    Vp2=np.reshape(Vp_arr,[ny,nz,nx])
    print('np.max(np.abs(Vp))='+str(np.max(np.abs(Vp2))))
    
    write_model_Vsp(rho,Vs2,Vp2)
    #time.sleep(20)
    os.system('mv %s' %'vs3dfile1.s Dump/vs3dfile_iter_'+str(it)+'.s')
    os.system('mv %s' %'vp3dfile1.p Dump/vp3dfile_iter_'+str(it)+'.p')
    os.system('mv %s' %'rho3dfile1.d Dump/rho3dfile_iter_'+str(it)+'.d')
    
    print('Finish Update iteration')
    #####################################
    #time.sleep(20)


    