#!/usr/bin/env python

from obspy.core import Stream, read, UTCDateTime
import glob
import operator
import numpy as np
from matplotlib.mlab import csd
from obspy.signal.invsim import paz_to_freq_resp
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import matplotlib.pyplot as plt
from math import pi
import math, copy
import matplotlib as mpl

mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=20)

# Here is the RS4D Nominal Response from OSOP

pazRS4D = {'gain': 6.56435, 'zeros': [-5779.15, 0, 0, 0], 'poles': [-163.344 + 102.457j,  
        -163.344 - 102.457j, -3.60818, -1.41864 + 0.410935j, -1.41864 - 0.410935j],
        'sensitivity': (335815000)}
        
pazRS4Dacc =  {'gain': 0.0141636, 'zeros': [657.269 + 1198.6j, 657.269 - 1198.6j, 0], 'poles': [-125.943 + 101.852j,  
        -125.943 - 101.852j, -0.0000623592],'sensitivity': (386825)}        
        
pazSTS2 = {'gain': 3.537344989*10**17,'zeros': [0., 0., -15.15, -176.6, -463.1 + 430.5j ,-463.1 - 430.5j], 'poles': [-0.037 - 0.037j,  
            -0.037 + 0.037j, -15.64, -97.34 - 400.7j, -97.34 + 400.7j, -374.8, -520.3, -10530 - 10050j, -10530 + 10050j, -13300, -255.097],
            'sensitivity': 20000.*0.9867629*((2**26)/40) ,'sensitivity': (1495)*((2.**24)/40.)*0.9867629}

pazTC = {'zeros': [0.j, 0.j, -392. + 0.j, -1960. + 0.j, -1490. + 1740.j, -1490. -1740.j],'poles': [-0.03691 + 0.03702j, -0.03691 - 0.03702j, -343. + 0.j, -370. + 467.j, -370. -467.j,
        -836. + 1522.j, -836. -1522.j, -4900. + 4700.j, -4900. - 4700.j, -6900. + 0.j, -15000. + 0.j],
        'gain': 4.344928*10**17, 'sensitivity': 751.4*0.9909*((2.**24)/40.)} 


# Start of Actual Script 


# set start and end times around the event
debug = True
day = 115
year = 2018
stime = UTCDateTime('2018-115T09:54:21.4')
etime = stime + 3.0


# load in the stations


stas = ['ASL9', 'FBA2', 'R141E']
locs = ['10', '10', '00']
chans = ['EHZ', 'EH0', 'EHZ']
sensors = ['STS2', 'Compact', 'RS4D']
net = ['GS', 'XX', 'AM']
color_list = ['k', 'green', 'm']
line_widths = ['4','2','2']



stF = Stream()


for slcs in zip(stas, locs, chans, sensors,net):
    
    files2 = glob.glob('/tr1/telemetry_days/' + slcs[4] + '_' + slcs[0] + '/' + str(year) + '/' + str(year) + '_' +  str(day).zfill(3) + '/' + slcs[1] + '_' + slcs[2] + '*')
    st = reduce(operator.add, map(read,files2))
    
    print(st)
    
    # You need to be careful as the clock gitter will make the trim function flaky
    
    
    if slcs[3] == 'STS2':
       paz = pazSTS2
       st.decimate(2,no_filter=True)
    elif slcs[3] == 'Compact':
        paz = pazTC
        st.decimate(2,no_filter=True)
    elif slcs[3] == 'RS4D':
        paz = pazRS4D    
        
    st.trim(stime-1,etime+1)
    delta = st[0].stats.delta
    if debug:
        print(delta)
    
        
    st.detrend('constant')
    st.taper(max_percentage=0.05, type='cosine')
    for tr in st:
        tr.simulate(paz_remove=paz,water_level=200)
        #tr.filter("bandpass",freqmin = 1.,freqmax= 20., corners=4.)
    
    stF += st
      
stF2 = stF.copy()           
            
 # Make the Figure 
 
 
idx = 0 
fig = plt.figure(1)
for tr in stF:
    tr.trim(stime, etime)
    t=np.arange(0, len(tr.data))/(tr.stats.sampling_rate)
    plt.plot(t,tr.data,linewidth=line_widths[idx],c=color_list[idx])
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (m/s)')
    plt.legend(['STS-2', 'Trillium Compact', 'RS4D geophone']) 
    idx=idx+1



# correct the 30 ms second timing issue 
count = 1
idx = 0
fig = plt.figure(2)
for tr in stF2:
    #tr.filter("bandpass",freqmin = 1.,freqmax= 20., corners=4.)
    tr.trim(stime, etime)
    if count == 3:
        t=np.arange(0, len(tr.data))/(tr.stats.sampling_rate)+0.00
    else:
        t=np.arange(0, len(tr.data))/(tr.stats.sampling_rate)
    plt.plot(t,tr.data, '.-', markersize=15,linewidth=line_widths[idx],c=color_list[idx])
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (m/s)')
    plt.legend(['STS-2', 'Trillium Compact', 'RS4D geophone']) 
    count = count+1 
    idx = idx+1
    




  
    
             
    
plt.show() 
