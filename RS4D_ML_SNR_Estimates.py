#!/usr/bin/env python
from obspy.core import read, UTCDateTime
from obspy.clients.fdsn import Client
import sys
import matplotlib.pyplot as plt
from obspy.signal.invsim import paz_to_freq_resp
from matplotlib.mlab import csd
import numpy as np
from obspy.io.xseed import Parser
from obspy.geodetics.base import gps2dist_azimuth



def rms(x):
    return np.sqrt(x.dot(x)/x.size)


def getmag(net, sta, sp, eve):
    # Function to get magnitudes
    eventtime = eve.origins[0].time
    
    print(sta.code)
    
    try:
    #if True:
        st= clientSH.get_waveforms('AM', sta.code , '00', 'SHZ', eventtime-120., eventtime+60., attach_response=True)
        print(st)
        
    except:
        try:
        #if True:
            st= clientSH.get_waveforms('AM', sta.code , '00', 'EHZ', eventtime-120., eventtime+60., attach_response=True)
            print(st)   
        
        except:
            print('No data for ' + sta.code + ' ' + str(eventtime))
            return 0.
        
        

    st.remove_response(output="DISP")
    for tr in st:
        #paz = sp.get_paz(tr.id, eventtime)
        #paz['zeros'].append(0.+0.j)
        tr.simulate(paz_simulate=paz_wa, water_level=100)
    print(st)
    #st = sp.rotate_to_zne(st)
    
    Noise = st.copy()
     
    # Compute Pre and post Event RMS values and get SNR
    
    #Noise = st.trim(eventtime-300, eventtime)
    
    st.trim(eventtime-5., eventtime+55.)
    Noise.trim(eventtime-120., eventtime-60.)
    
    st.plot()
    
    amp = []
    RMS_Noise = []
    RMS_Signal = []
    SNR = []
    
    
    
    for comp in ['Z']:
        amp.append(max(abs(st.select(component=comp)[0].data)))
        RMS_Noise.append(rms(Noise.select(component=comp)[0].data))
        RMS_Signal.append(rms(st.select(component=comp)[0].data))
        SNR.append(max(RMS_Signal) / max(RMS_Noise))
    
    
    # These are all lists -> convert everything to freaking numbers.     
    RMS_Noise = max(RMS_Noise)
    RMS_Signal = max(RMS_Signal)
    SNR = max(SNR)
        
        
    amp = max(amp)
    
    
    epi_dist, az, baz = gps2dist_azimuth(sta.latitude, sta.longitude,eve.origins[0].latitude,eve.origins[0].longitude)
    epi_dist = epi_dist / 1000
    
    print('Here is the amp:' + str(amp))
    a = 1.11
    b = 0.00189
    c = 2.09
    ml = np.log10(amp * 1E9) + a*np.log10(epi_dist) + b * epi_dist - c
    
    return ml, epi_dist, RMS_Noise, RMS_Signal, SNR


# Define the Response of a Wood-Anderson Sensor to Simulate what it would record
paz_wa = {'sensitivity': 1, 'zeros': [0 + 0j,0 + 0j], 'gain': 1.0028,
        'poles': [-5.49779 - 5.60886j, -5.49779 + 5.60886j]}

stime = UTCDateTime('2017-001T00:00:00.0')
etime = UTCDateTime('2018-091T00:00:00.0')
net = 'AM'

client = Client("IRIS")
clientSH = Client(base_url='https://fdsnws.raspberryshakedata.com/')

stas = clientSH.get_stations(maxlatitude=37.0, minlatitude= 35.0, minlongitude=-100, maxlongitude=-94.5, channel = "*HZ")
stas = stas[0]
for sta in stas:
    print(sta.code)
    print(sta)
    
    
compact_Ml = []
Neic_M = []

f=open('AM_SNR_Full2.txt','w')
f.write('Station, Calculated Magnitude, NEIC Magnitude, Event Distance, Event Lat, Event Lon,  RMS_Signal, RMS_Noise, SNR  \n')  
    
    
for sta in stas:
    try:

        coords = sta.longitude
    except:
        continue

    stalon = sta.longitude
    stalat = sta.latitude
    try:
    # For each station get the events within 2 degrees
        cat = client.get_events(starttime=stime, minmagnitude=2., latitude=stalat, 
                        longitude=stalon, maxradius=2.,endtime=etime)
    except:
        print('Network Down') 
        continue                   
    for eve in cat:
        for magt in eve.magnitudes:
            #print(magt)
            if magt.magnitude_type == 'ML':
                #print('In Loop')
                try:
                #if True:
                    ml,epi_dist,RMS_Noise,RMS_Signal,SNR  = getmag(net, sta, clientSH, eve)
                    print('ML=' + str(ml) + ' M=' + str(eve.magnitudes[0].mag))
                    compact_Ml.append(ml)
                    #print(compact_Ml)
                    Neic_M.append(eve.magnitudes[0].mag)

                    f.write(sta.code + ',' + str(ml) + ',' +
                    str(eve.magnitudes[0].mag) + ',' + str(epi_dist) + ',' +  
                    str(eve.origins[0].latitude) + ',' + str(eve.origins[0].longitude) + ',' + str(RMS_Noise) + ',' + str(RMS_Signal) + ',' + str(SNR)  + str( '\n'))
                except:
                    print('Could not compute ml')

f.close()

Neic_M = np.asarray(Neic_M)
compact_Ml = np.asarray(compact_Ml)

print(type(Neic_M))
Diff_Mag = compact_Ml - Neic_M



print(Diff_Mag)

fig = plt.figure(1)
plt.hist(Diff_Mag,300,facecolor='blue')
plt.ylabel('Counts')
plt.xlabel('Magnitude Misfit (GS-NEIC)')

plt.show()
            
            
        #except:
            #print('Could not compute ml')
