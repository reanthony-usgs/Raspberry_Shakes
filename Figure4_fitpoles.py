#!/usr/bin/env python
from obspy.core import read, UTCDateTime
import numpy as np
from scipy.optimize import fmin
import matplotlib.pyplot as plt
import math

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

paz = {'zeros': [0.j, 0.j, -392. + 0.j, -1960. + 0.j, -1490. + 1740.j, -1490. -1740.j],
    'poles': [-0.03691 + 0.03702j, -0.03691 - 0.03702j, -343. + 0.j, -370. + 467.j, -370. -467.j,
        -836. + 1522.j, -836. -1522.j, -4900. + 4700.j, -4900. - 4700.j, -6900. + 0.j, -15000. + 0.j],
        'gain': 4.344928*10**17, 'sensitivity': 754.3*2.**24/40.}
        


def corn_freq_2_paz(fc, damp=0.707):
    if damp <= 1.:
        poles = [-(damp + math.sqrt(1 - damp ** 2) * 1j) * 2 * np.pi * fc,
                 -(damp - math.sqrt(1 - damp ** 2) * 1j) * 2 * np.pi * fc]
    else:
        poles = [-(damp + math.sqrt(damp ** 2 - 1) * 1j) * 2 * np.pi * fc,
                 -(damp - math.sqrt(damp ** 2 - 1) * 1j) * 2 * np.pi * fc]
    return {'poles': poles, 'zeros': [0j, 0j], 'gain': 1, 'sensitivity': 1.0}


def mkresp(zp):
    f= zp[0]
    h = zp[1]
    #print(str(h) + ' '  +str(f))
    paz2 = corn_freq_2_paz(f,h)
    paz2['sensitivity']= zp[2]

    return paz2    
    
    

    
        

st1 = read('00_HH0.512.seed')
#print(st1)

st2 = read('00_EHZ.512.seed')
#print(st2)
stime = UTCDateTime('2018-156T20:15:34.1')
etime = stime + 60.


st1.trim(stime,etime)
st2.trim(stime,etime)
#st1.decimate(2)
#st1.decimate(2)
st1.simulate(paz_remove=paz)



def findresi(par):
    stTemp = st2.copy()
    paz2 = mkresp(par)
    stTemp.simulate(paz_remove=paz2)
    res = np.sum((stTemp[0].data-st1[0].data)**2)

    return res

#st2[0].data = -1*st2[0].data

st2.taper(.05)
st1.taper(.05)
st2.filter('bandpass',freqmin=.3, freqmax=2)
st1.filter('bandpass',freqmin=.3, freqmax=2)





norm= st2[0].std()/st1[0].std()

print('Here is the normalization: ' + str(norm))

#zp2 =[-154.6, 307, -37.5277, 0., 11.20, 0., 0.4448, 0., 301.4, 217.4, -33.81, 0., 10.9, 0., -2.564, 1.563]
#zp3 =[-542.8, 957.5, -167.7, 0., 10.72, 0., 0.4159, 0., -16.8, -63.3, -8.1, 0., -57.6, 0., 3.59, -6.28]
#zp4 =[-1568.3, 1782.5, -76.5, 0., 14.6, 0., 0.58, 0., -33.4, 5.652, -7.6434, 0., -225.41, 0.,  11.44, -23.144]

zp4=[.5, .707, 335815000.]
#print(mkresp(zp4))
bg = fmin(findresi,zp4, xtol=10**-6,ftol=10**-3,disp=False, maxiter=2000.)
#bg = zp4
print('Here is best fit: ' + str(bg) + ' for ' + st2[0].id)
paz2 = mkresp(bg)
#paz2 = mkresp(zp2)
#print(paz2)
st3 = st2.copy()
st3.simulate(paz_remove=mkresp(zp4))
st2.simulate(paz_remove=paz2)
t= np.arange(len(st1[0].data))/50.
fig = plt.figure(1,figsize=(12,12))
plt.subplots_adjust(hspace=0.001)
ax = plt.subplot(2,1,1)
ax.text(-0.07, 1.05, '(a)', transform=ax.transAxes,
      fontsize=18, fontweight='bold', va='top', ha='right')
plt.plot(t, st3[0].data, label='Nominal RS-4D Geophone', color='k')
plt.plot(t, st1[0].data, label='Reference Trillium Compact', color ='green' )
plt.plot(t, st2[0].data, label='Best-Fit RS-4D Geophone', color='m')
plt.legend(loc=1)
plt.xticks([])
plt.ylabel('Velocity (m/s)')
plt.ylim((-0.025, 0.025))
plt.xlim((min(t), max(t)))

    
ax=plt.subplot(2,1,2)
plt.plot(t, st3[0].data-st1[0].data, label='Nominal RS-4D Geophone', color='k')
#plt.plot(t, st1[0].data, label='Reference Trillium Compact', color ='green' )
plt.plot(t, st2[0].data- st1[0].data, label='Best-Fit RS-4D Geophone', color='m')
ax.text(-0.07, 1.05, '(b)', transform=ax.transAxes,
      fontsize=18, fontweight='bold', va='top', ha='right')
plt.ylabel('Residual Velocity (m/s)')
plt.xlabel('Time (s)')
plt.xlim((min(t), max(t)))
plt.ylim((-0.025, 0.025))
plt.legend(loc=1)
#plt.show()
plt.savefig('Result.eps', format='eps', dpi=400)
plt.savefig('Result.pdf', format='pdf', dpi=400)
