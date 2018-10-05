#!/usr/bin/env python

from obspy.core import read, UTCDateTime
import glob
import operator
import numpy as np
from matplotlib.mlab import csd
from obspy.signal.invsim import paz_to_freq_resp
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import matplotlib.pyplot as plt
from math import pi
import math, copy

#import matplotlib as mpl
#mpl.rc('font',family='serif')
#mpl.rc('font',serif='Times') 
#mpl.rc('text', usetex=True)
#mpl.rc('font',size=18)

# Compute Cross power
def cp(tr1, tr2, lenfft, lenol, delta):
    sr = 1./delta
    cpval, fre = csd(tr1.data, tr2.data,
                     NFFT=lenfft, Fs=sr, noverlap=lenol,
                     scale_by_freq=True)
    fre = fre[1:]
    cpval = cpval[1:]
    return cpval, fre

# Compute the response
def computeresp(resp, delta, lenfft):
    respval = paz_to_freq_resp(resp['poles'],
                               resp['zeros'],
                               resp['sensitivity']*resp['gain'],
                               t_samp=delta,
                               nfft=lenfft, freq=False)
    respval = np.absolute(respval*np.conjugate(respval))
    respval = respval[1:]
    return respval



    
def smoothF (frequency,power,samplingRate,octaveWindowWidth,octaveWindowShift,minFrequency,xStart):
   #
   ##########################################################################
   # smooth the spectra in the frequency domain
   #
   # To make the output independent of sample rate, the multiples are always with 
   # respect to XStart defined by user and not the Nyquist frequency
   #
   # HISTORY:
   #    2014-02-07 Manoch: created
   # 
   #
   ##########################################################################
   # 
   F     = [] # <-- frequency container
   P     = [] # <-- power container

   #
   # sort on frequency
   #
   frequency,power = (list(t) for t in zip(*sorted(zip(frequency,power))))

   halfWindow    = float(octaveWindowWidth / 2.0)
   windowShift   = octaveWindowShift

   #
   # maximum frequency at the Nyquist
   #
   maxFrequency = float(samplingRate) / 2.0

   #
   # xStart is set by user so that output remains independent of sample rate, the multiples are always with 
   # respect to xStart
   #
   Fc = xStart

   #
   # do the lower frequencies (<= xStart & >= minFrequency)
   #
   shift = math.pow(2.0,windowShift)
   while Fc >= minFrequency:

      thisBin = getBin(frequency,power,Fc,halfWindow)

      #
      # bin should not be empty
      #
      if (len(thisBin) > 0):
         P.append(np.mean(thisBin))
         F.append(Fc)
      else:
         P.append(float('NAN'))
         F.append(Fc)

      #
      # move the center frequency to the left by half of the windowWidth
      #
      Fc /= shift

   #
   # do higher frequencies (> xStart & <= Nyquist)
   #
   Fc = xStart * shift

   while Fc <= maxFrequency:

      thisBin = getBin(frequency,power,Fc,halfWindow)

      #
      # bin should not be empty
      #
      if (len(thisBin) > 0):
         P.append(np.mean(thisBin))
         F.append(Fc)
      else:
         P.append(float('NAN'))
         F.append(Fc)

      #
      # move the center frequency to the right by half of the windowWidth
      #
      Fc *= shift

   #
   # sort on frequency and return
   #
   F,P = (list(t) for t in zip(*sorted(zip(F,P))))
   return (F,P)    
    
def getBin(X,Y,Xc,octaveHalfWindow):
   #
   ##########################################################################
   # gather Y values that fall within octaveHalfWindow  on either side of the
   # centeral X value 
   # 
   # HISTORY:
   #    2014-02-07 Manoch: created
   # 
   ##########################################################################
   # 
   thisBin = []
   shift = math.pow(2.0,octaveHalfWindow)

   #
   # the bin is octaveHalfWindow around Xc
   #
   X1 = Xc / shift
   X2 = Xc * shift
   Xs = min(X1,X2)
   Xe = max(X1,X2)

   #
   # gather the values that fall within the range >=Xs and <= Xe
   #
   for i in xrange(0,len(X)):
      if X[i] >= Xs and X[i] <= Xe:
         thisBin.append(float(Y[i]))
   return thisBin   



   
#Compute self noise and PSD in dB
 
def compNoise(idx, length, overlap, delta, instresp):
    cpFix = lambda i1, i2: cp(st[i1], st[i2], length,
                              overlap, delta)

    # We could do these as permutations but that gets confusing
    # Instead we will just hard code the indices
    pp, f = cpFix(idx, idx)
    if idx == 0:
        noisetemp = pp - \
                    cpFix(1, 0)[0]*cpFix(0, 2)[0]/cpFix(1, 2)[0]
    elif idx == 1:
        noisetemp = pp - \
                    cpFix(2, 1)[0]*cpFix(1, 0)[0]/cpFix(2, 0)[0]
    elif idx == 2:
        noisetemp = pp - \
                    cpFix(1, 2)[0]*cpFix(2, 0)[0]/cpFix(1, 0)[0]
    else:
        print 'Bad index, crash landing.'
        sys.exit()
    # Convert to acceleration
    noisetemp *= (2.*np.pi*f)**2
    pp *= (2.*np.pi*f)**2
    # Remove the response
    noisetemp *= 1./instresp
    pp *= 1./instresp
    # Convert to dB
    noisetemp = 10.*np.log10(np.absolute(noisetemp))
    pp = 10.*np.log10(np.absolute(pp))
    return noisetemp, pp, f    
        


# Start of Actual Script 


debug = True
day = 72
year = 2018
stime = UTCDateTime('2018-072T03:00:00.0')
etime = stime + 60.*60.*6

#Here is my estimated response
#paz= {'poles': [(149.67189745697971+24.351503981794483j), (149.67189745697971-24.351503981794483j), -8.1574647449752771, -517.98353468592916, 
#    (-18.688964026369025-96.597388009454505j)], 'sensitivity': 495550.0, 
#    'zeros': [(-15324.27744543273+9331.9159286808681j), (-15324.27744543273-9331.9159286808681j), 
#    -653.22466280016897, 3.5520283683062983, 0.13662296274942703], 'gain': 1.0}

# Here is the RS4D Nominal Response from OSOP

pazRS4D = {'gain': 6.56435, 'zeros': [-5779, 0, 0, 0], 'poles': [-163.344 + 102.457j,  
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

files = glob.glob('/tr1/telemetry_days/AM_RB7A0/' + str(stime.year) + '/*'  +  str(stime.julday).zfill(3) + '/*EH*')
files += glob.glob('/tr1/telemetry_days/AM_R5D49/' + str(stime.year) + '/*'  +  str(stime.julday).zfill(3) + '/*EH*')
files += glob.glob('/tr1/telemetry_days/AM_R141E/' + str(stime.year) + '/*'  +  str(stime.julday).zfill(3) + '/*EH*')


st = reduce(operator.add, map(read,files))
# You need to be careful as the clock gitter will make the trim function flaky
st.trim(starttime=stime, endtime=etime)
if debug:
    print(st)

# Now we have our data why not compute some spectra and self noise?
length =2**18
overlap=2**14
delta = st[0].stats.delta
instresp = computeresp(pazRS4D, delta, length)

n1,p1,fS = compNoise(0, length, overlap, delta, instresp)
n2,p2,fS = compNoise(1, length, overlap, delta, instresp)
n3,p3,fS = compNoise(2, length, overlap, delta, instresp)

Noise_Mean = (n1+n2+n3)/3.
Power_Mean = (p1+p2+p3)/3.


Smoothing_Window = 1/64.;

fSm, S_NM = smoothF(fS,Noise_Mean,1/delta,Smoothing_Window,Smoothing_Window*2,0.01,50.)
fSm, S_PM = smoothF(fS,Power_Mean,1/delta,Smoothing_Window,Smoothing_Window*2,0.01,50.)

#fSm, S_n1 = smoothF(fS,n1,1/delta,Smoothing_Window,Smoothing_Window*2,0.01,50.)
#fSm, S_n2 = smoothF(fS,n2,1/delta,Smoothing_Window,Smoothing_Window*2,0.01,50.)
#fSm, S_n3 = smoothF(fS,n3,1/delta,Smoothing_Window,Smoothing_Window*2,0.01,50.)


fSm = np.array(fSm)
S_NM = np.array(S_NM)
S_PM = np.array(S_PM)

#S_n1 = np.array(S_n1)
#S_n2 = np.array(S_n2)
#S_n3 = np.array(S_n3)



NLNMper,NLNMpower = get_nlnm()
NHNMper,NHNMpower = get_nhnm()

# Now load in the STS2, Trillium Compact, and Shake Accelerometer Data and Get Plot the Acceleration Spectra

Smoothing_Window = 1/64.;

stas = ['ASL9', 'FBA2', 'R5D49']
locs = ['10', '10', '00']
chans = ['EHZ', 'EH0', 'ENZ']
sensors = ['STS2', 'Compact', 'RS4D-Accel']
net = ['GS', 'XX', 'AM']



for slcs in zip(stas, locs, chans, sensors,net):
    
    files2 = glob.glob('/tr1/telemetry_days/' + slcs[4] + '_' + slcs[0] + '/' + str(year) + '/' + str(year) + '_' +  str(day).zfill(3) + '/' + slcs[1] + '_' + slcs[2] + '*')
    st = reduce(operator.add, map(read,files2))
    
    # You need to be careful as the clock gitter will make the trim function flaky
    st.trim(starttime=stime, endtime=etime)
    if debug:
        print(st)
    if slcs[3] == 'STS2':
       paz = pazSTS2
       st.decimate(2)
    elif slcs[3] == 'Compact':
        paz = pazTC
        st.decimate(2)
    elif slcs[3] == 'RS4D-Accel':
        paz = pazRS4Dacc    
        
    st.trim(stime, etime)
    delta = st[0].stats.delta
    if debug:
        print(delta)
        
    instresp = computeresp(paz,delta,length)
    
    for tr in st:
        (p,f) = cp(tr,tr,length,overlap,delta)
        p = p/instresp
        fSm2, pSm = smoothF(f,p,1/delta,Smoothing_Window,Smoothing_Window*2,0.01,50.)
        fSm2 = np.array(fSm2)
        pSm = np.array(pSm)
        
        
        if slcs[3] == 'STS2':
            psd = 10.*np.log10(((2*pi*fSm2)**2)*pSm)
            psd = psd.real
            PSDSTS2 = psd
            fSTS2 =fSm2
            
        elif slcs[3] == 'Compact':
            psd = 10.*np.log10(((2*pi*fSm2)**2)*pSm)
            psd = psd.real
            PSDTC = psd
            fTC = fSm2
            
            
        elif slcs[3] == 'RS4D-Accel':
            psd = 10.*np.log10(((2*pi*fSm2)**0)*pSm)
            psd = psd.real
            PSDRSA = psd
            fRSA = fSm2    
            
            
           
            
NLNMper,NLNMpower = get_nlnm()
NHNMper,NHNMpower = get_nhnm()


# Write all the data out to textfiles

import csv
f=open('ShakesNoise.txt','w')

writer = csv.writer(f, delimiter='\t')
writer.writerows(zip(fSm,S_NM,S_PM))

f.close()

f=open('PSD_Compact.txt','w')

writer = csv.writer(f, delimiter='\t')
writer.writerows(zip(fTC,PSDTC))

f.close()


f=open('PSD_STS2.txt','w')

writer = csv.writer(f, delimiter='\t')
writer.writerows(zip(fSTS2,PSDSTS2))

f.close()


f=open('PSD_MEMS.txt','w')

writer = csv.writer(f, delimiter='\t')
writer.writerows(zip(fRSA,PSDRSA))

f.close()






fig =plt.figure(1,figsize=(12,12))
#plt.semilogx(fSm,S_n1,'grey', linewidth=2.,label= 'Raspberry Shake 4D (Geophone) Self-Noise')
#plt.semilogx(fSm,S_n2,'dimgrey', linewidth=2.)
#plt.semilogx(fSm,S_n3,'slategrey', linewidth=2.)
plt.semilogx(fSm,S_NM,'grey', linewidth=2.,label= 'RS4D (Geophone) Self-Noise')
plt.semilogx(fSm,S_PM, 'r', linewidth=2.,label = 'RS4D (Geophone) PSD')
plt.semilogx(fRSA,PSDRSA, 'm', linewidth=2.,label = 'RS4D (Accelerometer) PSD/Self-Noise')
plt.semilogx(fTC,PSDTC, 'g', linewidth=2.,label = 'Trillium Compact PSD')
plt.semilogx(fSTS2,PSDSTS2, 'b', linewidth=2.,label = 'STS-2 PSD')
plt.semilogx(1./NLNMper, NLNMpower, linewidth=4., color='k')
plt.semilogx(1./NHNMper, NHNMpower, linewidth=4., color='k',label='NLNM/NHNM')
plt.xlim((50.,.0333))
#plt.ylim((-180., -80.))
plt.ylabel('Power (dB rel. 1 $(m/s^2)^2/Hz)$')
plt.xlabel('Frequency (Hz)')
#plt.title('Self-Noise Estimates for ' + str(stime.year) + ' ' + str(stime.julday).zfill(3) + str(stime.hour).zfill(2) + ':' + str(stime.minute).zfill(2) + ' Duration 2 Hour')
plt.legend()
plt.show()   





    
            
            
           
        
    














         
        
        
         


