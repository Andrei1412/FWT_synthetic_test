#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 14:19:45 2018

@author: user
"""

import numpy as np
#import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter
from scipy import signal
from qcore import timeseries
import os

def readGP_2(loc, fname):
    """
    Convinience function for reading files in the Graves and Pitarka format
    """
    with open("/".join([loc, fname]), 'r') as f:
        lines = f.readlines()
    
    data = []

    for line in lines[2:]:
        data.append([float(val) for val in line.split()])

    data=np.concatenate(data) 
    
    line1=lines[1].split()
    num_pts=float(line1[0])
    dt=float(line1[1])
    shift=float(line1[4])

    return data, num_pts, dt, shift

def computeFourier(accTimeSeries, dt, duration):
    #computes fourier spectra for acceleration time series
    # TODO: compute original number of points (setting to default for now) change to npts = len(accTimeSeries)
    npts = len(accTimeSeries)
    npts_FFT = int(np.ceil(duration)/dt)
    
    # compute number of points for efficient FFT
    ft_len = int(2.0 ** np.ceil(np.log(npts_FFT) / np.log(2.0)))
    if npts > ft_len:
        accTimeSeries = accTimeSeries[:ft_len]
        npts = len(accTimeSeries)
    
    # Apply hanning taper to last 5% of motion
    ntap = int(npts * 0.05) 
    accTimeSeries[npts - ntap:] *= np.hanning(ntap * 2 + 1)[ntap + 1:]
    
    # increase time series length with zeroes for FFT
    accForFFT = np.pad(accTimeSeries, (0, ft_len - len(accTimeSeries)), 'constant', constant_values=(0,0))
    ft = np.fft.rfft(accForFFT)
    
    # compute frequencies at which fourier amplitudes are computed
    ft_freq = np.arange(0, ft_len / 2 + 1) * ( 1.0 / (ft_len * dt))
    return ft, ft_freq

def butter_bandpass(lowcut, highcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

######################################
def normpdf_python(x, mu, sigma):
   return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-1*(x-mu)**2/2*sigma**2)

#def Calc_Displacement(Vel,Vel_O,t):
#    nt=len(t)
#    dt=t[1]-t[0]
#    inc_z=0
#    inc_zo=0
#    Dis_S=np.zeros(Vel.shape)
#    Dis_O=np.zeros(Vel.shape)
#    for m in np.arange(1,nt):
#        inc_z+= Vel[m]
#        Dis_S[m]=inc_z*dt
#        inc_zo+= Vel_O[m]
#        Dis_O[m]=inc_zo*dt        
#    return Dis_S, Dis_O
    
##############################################
def source_adj_1(stat_data_S,stat_data_O,num_pts,dt,b, a):

#    t = np.arange(num_pts)*dt
#    Dis_S, Dis_O = Calc_Displacement(stat_data_S,stat_data_O,t)
    Dis_S = np.cumsum(stat_data_S)*dt
    Dis_O = np.cumsum(stat_data_O)*dt    
    
    Dis_S  = np.multiply(signal.tukey(int(num_pts),0.1),Dis_S)    
    Dis_O  = np.multiply(signal.tukey(int(num_pts),0.1),Dis_O)

    Dis_S = signal.filtfilt(b, a, Dis_S)
    Dis_O = signal.filtfilt(b, a, Dis_O)
    
    Jp = Dis_S - Dis_O
    Source = np.flip(Jp, axis=0)

    return Source

def write_adj_source(s1,v1,mainfolder,mainfolder_source,source):
    
    with open("/".join([mainfolder, s1]), 'r') as f:
        lines = f.readlines()
    tline1=   lines[0]     
    tline2=   lines[1]
    
    filename1=mainfolder_source+v1
    print(filename1)
    fid = open(filename1,'w')
    fid.write("%s" %(tline1))
    fid.write("%s" %(tline2))
    lt=len(source)
    count=0
    while (count+1)*6<lt:
        fid.write("%10f%10f%10f%10f%10f%10f%s" %(source[count*6],source[count*6+1],source[count*6+2],source[count*6+3],source[count*6+4],source[count*6+5],'\n'))
        count+=1
    ii=lt-count*6   
    i=0
    while (i<ii):
        i+=1
        fid.write("%10f%s" %(source[lt-ii+i-1],'\n'))
    fid.close()
    return

def write_adj_source_ts(s1,v1,mainfolder,mainfolder_source,source,dt):
    #filename1=mainfolder_source+v1
    vs1=v1.split('.')
    timeseries.seis2txt(source,dt,mainfolder_source,vs1[0],vs1[1])
    return	


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
    
###################################################
#statnames = ['AAA','BBB','CCC','DDD','EEE','FFF','GGG','HHH','III','JJJ','KKK','LLL','MMM','NNN','PPP','QQQ','RRR','SSS','TTT','UUU','VVV','XXX','YYY','ZZZ']
#statnames = ['AAA','BBB','CCC','DDD']
#statnames = ['A1' ,'A2' ,'A3' ,'A4' ,'B1' ,'B2' ,'B3' ,'B4' ,'C1' ,'C2' ,'C3' ,'C4' ,'D1' ,'D2' ,'D3' ,'D4']
statnames = ['A1' ,'A2' ,'A3' ,'A4' ,'A5' ,'A6' , 'A7', 'B1' ,'B2' ,'B3' ,'B4' ,'B5' ,'B6' , 'B7','C1' ,'C2' ,'C3' ,'C4' ,'C5' ,'C6' ,'C7','D1' ,'D2' ,'D3' ,'D4' ,'D5' ,'D6' ,'D7','E1' ,'E2' ,'E3' ,'E4' ,'E5' ,'E6' ,'E7','F1' ,'F2' ,'F3' ,'F4' ,'F5' ,'F6','F7','G1' ,'G2' ,'G3' ,'G4' ,'G5' ,'G6','G7']

#statnames = ['A1' ,'A2' ,'A3' ,'A4' ,'A5', 'B1' ,'B2' ,'B3' ,'B4' ,'B5', 'C1' ,'C2' ,'C3' ,'C4','C5' ,'D1' ,'D2' ,'D3' ,'D4','D5','E1' ,'E2' ,'E3' ,'E4','E5']
GV=['.090','.000','.ver']
GV_ascii=['.x','.y','.z']

mainfolder='../../Vel_es/Vel_es_i/'
mainfolder_o='../../Vel_ob/Vel_ob_i/'
mainfolder_source='../../../AdjSims/V3.0.7-a2a_xyz/Adj-InputAscii/'
os.system('rm ../../../AdjSims/V3.0.7-a2a_xyz/Adj-InputAscii/*.*')

print(mainfolder_o)
_, num_pts, dt, shift  = readGP_2('../../Vel_ob/Vel_ob_i','A1.000')
t = np.arange(num_pts)*dt
############/nesi/nobackup/nesi00213/RunFolder/tdn27/rgraves/Adjoint/Syn_VMs/Kernels/#########################
fs = 1/dt
lowcut = 0.01
highcut = 0.2

fc = highcut  # Cut-off frequency of the filter
w = fc / (fs / 2) # Normalize the frequency
b, a = signal.butter(4, w, 'low')
#inv_data = butter_bandpass_filter(data, lowcut, highcut, fs, order=4)
source_file='../../../StatInfo/SOURCE.txt'
station_file='../../../StatInfo/STATION.txt'

nShot, S = read_source(source_file)
nRec, R = read_source(station_file)
################################
#distance[j-1]=((R[j-1,1]-S[1])**2+(R[j-1,2]-S[2])**2+(R[j-1,0]-S[0])**2)**(0.5)
fi1=open('iShot.dat','r')
ishot=np.int64(np.fromfile(fi1,dtype='int64'))
fi1.close()
print('ishot='+str(ishot))

for i,statname in enumerate(statnames):

    s0=statname+GV[0]
    v0=statname+GV_ascii[0]
#    
    s1=statname+GV[1]
    v1=statname+GV_ascii[1]
#
    s2=statname+GV[2]
    v2=statname+GV_ascii[2]    
#    
   
    #stat_data_0_S, num_pts, dt, shift  = readGP_2(mainfolder,s0)
    #stat_data_1_S, _, _, _  = readGP_2(mainfolder,s1)
    #stat_data_2_S, _, _, _  = readGP_2(mainfolder,s2)
    stat_data_0_S  = timeseries.read_ascii(mainfolder+s0)
    stat_data_0_S  = np.multiply(signal.tukey(int(num_pts),0.1),stat_data_0_S)
    stat_data_1_S  = timeseries.read_ascii(mainfolder+s1)
    stat_data_1_S  = np.multiply(signal.tukey(int(num_pts),0.1),stat_data_1_S)   
    stat_data_2_S  = timeseries.read_ascii(mainfolder+s2)
    stat_data_2_S  = np.multiply(signal.tukey(int(num_pts),0.1),stat_data_2_S)  
    #stat_data_0_O, num_pts, dt, shift1  = readGP_2(mainfolder_o,s0)
    #stat_data_1_O, _, _, _  = readGP_2(mainfolder_o,s1)
    #stat_data_2_O, _, _, _  = readGP_2(mainfolder_o,s2)    
    stat_data_0_O  = timeseries.read_ascii(mainfolder_o+s0)
    stat_data_0_O  = np.multiply(signal.tukey(int(num_pts),0.1),stat_data_0_O)
    stat_data_1_O  = timeseries.read_ascii(mainfolder_o+s1)
    stat_data_1_O  = np.multiply(signal.tukey(int(num_pts),0.1),stat_data_1_O)
    stat_data_2_O  = timeseries.read_ascii(mainfolder_o+s2)
    stat_data_2_O  = np.multiply(signal.tukey(int(num_pts),0.1),stat_data_2_O)   


    
    t = np.arange(num_pts)*dt
    nt=num_pts
        
    stat_data_0_Sf = stat_data_0_S
    stat_data_1_Sf = stat_data_1_S
    stat_data_2_Sf = stat_data_2_S
     
    stat_data_0_Of = stat_data_0_O
    stat_data_1_Of = stat_data_1_O
    stat_data_2_Of = stat_data_2_O
       
#    stat_data_0_Sf = signal.filtfilt(b, a, stat_data_0_S)
#    stat_data_1_Sf = signal.filtfilt(b, a, stat_data_1_S)
#    stat_data_2_Sf = signal.filtfilt(b, a, stat_data_2_S)
#    
#    stat_data_0_Of = signal.filtfilt(b, a, stat_data_0_O)
#    stat_data_1_Of = signal.filtfilt(b, a, stat_data_1_O)
#    stat_data_2_Of = signal.filtfilt(b, a, stat_data_2_O)    
    
    print('ireceiver='+str(i))
    distance=((R[i,1]-S[ishot-1,1])**2+(R[i,2]-S[ishot-1,2])**2+(R[i,0]-S[ishot-1,0])**2)**(0.5)
    
    if((distance<500) and (distance>50)):
        source_x=source_adj_1(stat_data_0_Sf,stat_data_0_Of,num_pts,dt,b, a)
        source_y=source_adj_1(stat_data_1_Sf,stat_data_1_Of,num_pts,dt,b, a)
#        source_x=np.zeros(stat_data_0_Sf.shape)
#        source_y=np.zeros(stat_data_0_Sf.shape)        
        source_z=source_adj_1(stat_data_2_Sf,stat_data_2_Of,num_pts,dt,b, a)   
    else:
        source_x=np.zeros(stat_data_0_Sf.shape)
        source_y=np.zeros(stat_data_0_Sf.shape)
        source_z=np.zeros(stat_data_0_Sf.shape)

    write_adj_source_ts(s0,v0,mainfolder,mainfolder_source,source_x,dt)
    write_adj_source_ts(s1,v1,mainfolder,mainfolder_source,source_y,dt)
    write_adj_source_ts(s2,v2,mainfolder,mainfolder_source,source_z,dt)    
    
#f_err = open('err_shot1.dat','w')
#Err.astype('float').tofile(f_err)     
    
    
    
    
