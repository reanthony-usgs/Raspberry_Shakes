#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as ml
from mpl_toolkits.basemap import Basemap, maskoceans
import matplotlib.pyplot as plt
jet = plt.cm.get_cmap('jet')
bwr = plt.cm.get_cmap('coolwarm')


#import cartopy
#import cartopy.io.shapereader as shpreader
#import cartopy.crs as ccrs


debug = False

#get the station information 

stas = []
lats = []
lons = []

with open('Shakes_NEIC_misfit.csv') as f:
    next(f)
    for line in f:
        line = line.split(',')
        stas.append(line[1])
        lats.append(float(line[2]))
        lons.append(float(line[3]))
f.close()


# get the Event information 

EQ_M = []
EQ_lats = [] 
EQ_lons = []

with open('T_Compacts_Z.txt') as f:
    next(f)
    for line in f:
        line = line.split(',')
        EQ_M.append(float(line[2]))
        EQ_lats.append(float(line[4]))
        EQ_lons.append(float(line[5]))
f.close()

with open('Event_Magnitudes_AM_Full.txt') as f:
    next(f)
    for line in f:
       for line in f:
        line = line.split(',')
        EQ_M.append(float(line[2]))
        EQ_lats.append(float(line[4]))
        EQ_lons.append(float(line[5]))
f.close()

print(len(EQ_M))
print(type(EQ_M))
print(EQ_M[0])


# Now we have all the EQ and station locations - make the map

print('Map parameters', min(EQ_lats)-0.25, max(EQ_lats)+0.25, min(EQ_lons)-0.25, max(EQ_lons)+1.)



# Now we make the plot
fig = plt.figure(1,figsize=(12,6)) 
m = Basemap(projection='merc',llcrnrlat=min(EQ_lats)-0.25,urcrnrlat=max(EQ_lats)+0.25,\
           llcrnrlon=min(EQ_lons)-0.25,urcrnrlon=max(EQ_lons)+1.,lat_ts=20,resolution='h',area_thresh=500.)


m.drawmapboundary(fill_color='c')
m.fillcontinents(color='beige', alpha=1.0, lake_color='w')

min_size_ = min(EQ_M) -1
max_size_ = max(EQ_M) +1
print(max_size_)
print(min_size_)
frac = [(0.01 + (_i - min_size_)) / (max_size_ - min_size_) for _i in EQ_M]
size_plot = [8*(_i * (max_size_ - min_size_)) ** 4 for _i in frac]



xG, yG = m(EQ_lons, EQ_lats)
sc = plt.scatter(xG, yG, s=size_plot, alpha=0.5, marker='o', edgecolors='white', zorder=2)



# plot the scale for EQ Magnitude 
mags2 = [2., 2.5, 3., 3.5, 4.]
frac2 = [(0.01 + (_i - min_size_)) / (max_size_ - min_size_) for _i in mags2]
size_plot2 = [8*(_i * (max_size_ - min_size_)) ** 4 for _i in frac2]


for mag, siz, f2 in zip(mags2,size_plot2, frac2):
    sc2 = m.scatter([] ,[], alpha=0.8, s=siz, c='#1f77b4', zorder=3,
    label='$M_L$' + str(mag))
    EQ_L = plt.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0.)

m.drawcoastlines(linewidth=2, zorder=3)
m.drawcountries(zorder=3)
m.drawstates(zorder=3)


# Plot the Seismic Stations
xGs, yGs = m(lons, lats)
sc = plt.scatter(xGs[0:17], yGs[0:17], alpha=1.0, s=150, marker='^',color='forestgreen', edgecolors='black', label='USGS Aftershock', zorder=4)
sc3 = plt.scatter(xGs[17:40], yGs[17:40], alpha=1.0, s=150, marker='^',color='m', edgecolors='black', label='Raspberry Shake',zorder=4)


leg = plt.legend(handles=[sc,sc3],fancybox=True, loc='best',fontsize=10)
leg.get_frame().set_alpha(0.93)




m.drawparallels(np.arange(35,40,1),linewidth=0.33,labels=[1,1,0,0],zorder=4)
m.drawmeridians(np.arange(-94,-104,-2), linewidth=0.33,labels=[0,0,0,1],zorder=4)

from matplotlib import pyplot
pyplot.gca().add_artist(EQ_L)

# Plot Major Cities (OCK,Tulsa,Enid)

C_Lats = [35.4676, 36.1540, 36.3956, 35.2226]
C_Lons = [-97.5164, -95.9928, -97.8784, -97.4395]

xCL, yCL = m(C_Lons, C_Lats)

cc = plt.scatter(xCL[0:4],yCL[0:4], marker='o',color='black',s=50,zorder=5)
#cc = plt.scatter(xCL[1],yCL[1], marker='o',color='black',s=50,zorder=5)
#cc = plt.scatter(xCL[2],yCL[2], marker='o',color='black',s=50,zorder=5)


fig.tight_layout()

plt.show()
