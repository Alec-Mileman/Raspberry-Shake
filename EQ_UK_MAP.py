#========================================================
#========================================================
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd
#========================================================
# This section opens up the 2 text files 'UK_EQ_STATIONS.txt' and 'UK_EQ.txt',
# Reading them into dataframes for use in the code

with open('UK_EQ_STATIONS.txt', 'r') as logs:
    logs.readline()
    columns = ['Network', 'Station', 'Latitude', 'Longitude', 'Elevation', 'Sitename', 'StartTime', 'EndTime']
    formatted_logs = pd.DataFrame([dict(zip(columns, line.split('|'))) for line in logs])

stations_lat = formatted_logs['Latitude'][3:]
stations_lons = formatted_logs['Longitude'][3:]

lons_stations = stations_lons.values.tolist()
lats_stations = stations_lat.values.tolist()

for i in range(len(stations_lons)):
    lons_stations[i] = float(lons_stations[i])
    lats_stations[i] = float(lats_stations[i])

#========================================================
with open('UK_EQ.txt', 'r') as logs_eq:
    logs_eq.readline()
    columns = ['DATE', 'LATITUDE', 'LONGITUDE', 'MAGNITUDE', 'DEPTH', 'TIME']
    formatted_logs_eq = pd.DataFrame([dict(zip(columns, line.split(','))) for line in logs_eq])

eq_lat_1 = formatted_logs_eq['LATITUDE']
eq_lons_1 = formatted_logs_eq['LONGITUDE']
eq_mag_1 = formatted_logs_eq['MAGNITUDE']

eq_lats = eq_lat_1.values.tolist()
eq_lons = eq_lons_1.values.tolist()
eq_mag = eq_mag_1.values.tolist()

for i in range(len(eq_lats)):
    eq_lats[i] = float(eq_lats[i])
    eq_lons[i] = float(eq_lons[i])

for i in range(len(eq_mag)):
    eq_mag[i] = float(eq_mag[i])

#========================================================
# This sets the marker size for the EQs, organising their ML into marker sizes

marker_size = []

for i in range(len(eq_mag)):
     if 2.0 <= eq_mag[i] < 2.5:
         marker_size.append(10)
         
     elif 2.5 <= eq_mag [i] < 3.0:
         marker_size.append(15)
         
     elif 3.0 <= eq_mag[i] < 3.5:
        marker_size.append(20)
        
#========================================================

fig = plt.figure(figsize=(20,20))
res = 4000 # This sets the resolution of the image. I found 4000 to be a good quality. Lower for the code to run faster.


# Setting up map image through basemap
m = Basemap(llcrnrlon = -10.5544, llcrnrlat = 49.35, # Area of domain for map plotting, this is for the UK
            urcrnrlon = 3, urcrnrlat = 60.905,
            resolution = 'c', # Set using letters, e.g. c is a crude drawing, f is a full detailed drawing
            epsg = 27700,
            projection = 'merc', # The projection style is what gives us a 2D view of the world for this
            lon_0 = -4.36,lat_0 = 54.7) # Setting the central point of the image
            
# Downloading an image of the UK to be fit within the map
m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = res, verbose= True)

# This sets the longitudes and latitudes of the stations, from the df created about.
# This is used to plot the station markers in the correct positions
x_1, y_1 = m(lons_stations, lats_stations)

# The same as above, but for EQ locations
x_eq, y_eq = m(eq_lons, eq_lats)

# Indiviudally plotting the EQ locations and their corresponding marker size
for x in range(len(marker_size)):
    m.plot(x_eq[x], y_eq[x],'o', color = 'r', markersize = marker_size[x])

# Plotting station locations
m.plot(x_1, y_1, '^', color = 'w', markersize = 7, label = 'RSPSHK Stations')

#A janky way of seting the legend markers
m.plot(x_eq[0], y_eq[0], 'o', color = 'r', markersize = 10, label = '2.0 - 2.5 ML')
m.plot(x_eq[0], y_eq[0], 'o', color = 'r', markersize = 15, label = '2.5 - 3.0 ML')
m.plot(x_eq[0], y_eq[0], 'o', color = 'r', markersize = 20, label = '3.0 - 3.5 ML')

# Setting the meridians
parallels = np.arange(0.,81,1.)
m.drawparallels(parallels,labels=[False,True,True,False], fontsize = 20)
meridians = np.arange(-351.,351.,1.)
m.drawmeridians(meridians,labels=[True,False,False,True], fontsize = 20)

#========================================================

plt.legend(loc = 'upper left', prop = {'size':25})
plt.show()
