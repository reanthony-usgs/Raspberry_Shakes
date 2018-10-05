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

def getstalist(sp, etime, net, debug=False):
    # Inputs - Decoded dataless using Parser function ,end-time, Network
    """ A function to get a station list. """
    stations = []
    for cursta in sp.stations:
# As we scan through blockettes we need to find blockettes 50
        for blkt in cursta:
            if blkt.id == 50:
# Pull the station info for blockette 50
                stacall = blkt.station_call_letters.strip()
                if debug:
                    print "Here is a station in the dataless: " + stacall
                if type(blkt.end_effective_date) is str:
                    curdoy = strftime("%j", gmtime())
                    curyear = strftime("%Y", gmtime())
                    curtime = UTCDateTime(curyear + "-" + curdoy +"T00:00:00.0")

                    if blkt.start_effective_date <= etime:
                        stations.append(blkt.station_call_letters.strip())
                elif blkt.start_effective_date <= etime and blkt.end_effective_date >= etime:
                    stations.append(blkt.station_call_letters.strip())
    return stations


def getmag(net, sta, sp, eve):
    # Function to get magnitudes - Input Seismic Network, Stations, ???, Events
    eventtime = eve.origins[0].time
    try:
        st = read('/msd/' + net + '_' + sta + '/' + str(eventtime.year) +
                  '/' + str(eventtime.julday).zfill(3) + '/*HH*')
    except:
        print('No data for ' + sta + ' ' + str(eventtime))
        return 0.

    st.trim(eventtime-120., eventtime+60.)
    
    #Noise = st

    for tr in st:
        paz = sp.get_paz(tr.id, eventtime)
        #print(paz)
        paz['zeros'].append(0.+0.j)
        tr.simulate(paz_remove=paz, paz_simulate=paz_wa, water_level=100)

    print(st)
    st = sp.rotate_to_zne(st)
    
    Noise = st.copy()
    
    # Compute Pre and post Event RMS values and get SNR
    
    #Noise = st.trim(eventtime-300, eventtime)
    
    st.trim(eventtime-5., eventtime+55.)
    Noise.trim(eventtime-120., eventtime-60.)
    
    
    amp = []
    RMS_Noise = []
    RMS_Signal = []
    SNR = []
    
    for comp in ['Z']:
        amp.append(max(abs(st.select(component=comp)[0].data)))
        RMS_Noise.append(rms(Noise.select(component=comp)[0].data))
        print(RMS_Noise)
        RMS_Signal.append(rms(st.select(component=comp)[0].data))
        SNR.append(max(RMS_Signal) / max(RMS_Noise))
    
    # These are all lists -> convert everything to freaking numbers.     
    RMS_Noise = max(RMS_Noise)
    RMS_Signal = max(RMS_Signal)
    SNR = max(SNR)
    
        
    amp = max(amp)
    print(SNR)

    # Get station coordinates
    coords = sp.get_coordinates(net + '.' + sta + '.00.HHZ', stime)
    # Map of events

    epi_dist, az, baz = gps2dist_azimuth(coords['latitude'],
    coords['longitude'],eve.origins[0].latitude,eve.origins[0].longitude)
    
    # Convert distance to EpiCenter to km
    epi_dist = epi_dist / 1000.

    print(epi_dist)

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


net = 'GS'

client = Client("IRIS")

# Obspy Function to read to dataless for an entire Network
sp = Parser('/APPS/metadata/SEED/' + net + '.dataless')

# Get oklahoma stations
stas = getstalist(sp, stime, net)
goodstas = []
for sta in stas:
    # Find only Oklahoma Aftershock Deployments
    if 'OK0' in sta:
        goodstas.append(sta)
stas = goodstas

compact_Ml = []
Neic_M = []

f=open('SNR_T_Compacts.txt','w')
f.write('Station, Calculated Magnitude, NEIC Magnitude, Event Distance, Event Lat, Event Lon, RMS_Signal, RMS_Noise, SNR \n')


for sta in stas:
    try:
        coords = sp.get_coordinates(net + '.' + sta + '.00.HHZ', stime)
    except:
        continue

    # Map of events
    stalon = coords['longitude']
    stalat = coords['latitude']
    # For each station get the events within 2 degrees
    
    try:
        cat = client.get_events(starttime=stime, minmagnitude=2., latitude=stalat,
                        longitude=stalon, maxradius=2.,endtime=etime)
    except:
        print('Connection Reset')                    
    for eve in cat:                    
        for magt in eve.magnitudes:
            #print(magt)
            if magt.magnitude_type == 'ML':
                #print('In Loop')
                try:
                #if True:                         
    
                    ml,epi_dist,RMS_Noise,RMS_Signal,SNR = getmag(net, sta, sp, eve)

                    print('ML=' + str(ml) + ' ' + eve.magnitudes[0].magnitude_type + '=' + str(eve.magnitudes[0].mag))
                    compact_Ml.append(ml)

                    Neic_M.append(eve.magnitudes[0].mag)

                    f.write(sta + ',' + str(ml) + ',' +
                    str(eve.magnitudes[0].mag) + ',' + str(epi_dist) + ',' +  
                    str(eve.origins[0].latitude) + ',' + str(eve.origins[0].longitude) + ',' + str(RMS_Noise) + ',' + str(RMS_Signal) + ',' + str(SNR) + str( '\n'))

                except:
                    print('Could not compute ml')

f.close()

Neic_M = np.asarray(Neic_M)
compact_Ml = np.asarray(compact_Ml)

print(type(Neic_M))
Diff_Mag = compact_Ml - Neic_M



print(Diff_Mag)

fig = plt.figure(1)
plt.hist(Diff_Mag,30,facecolor='blue')
plt.ylabel('Counts')
plt.xlabel('Magnitude Misfit (GS-NEIC)')

plt.show()

# Cat all of these things into an array to save

# np.savetxt('Cmpact_Magnitudes_2_deg_120_130.txt', Neic_M, fmt='%.3f')
