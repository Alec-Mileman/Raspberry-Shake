#============================================================
#============================================================
import matplotlib.pyplot as plt
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn import Client

para = {'axes.labelsize': 20, 'font.size': 20, 'legend.fontsize': 16,
'xtick.labelsize': 20,'ytick.labelsize': 20, 'figure.subplot.left': 0.12,
'figure.subplot.right': 0.98, 'figure.subplot.bottom': 0.11,
'figure.subplot.top': 0.97,'axes.linewidth': 1.0}

plt.rcParams.update(para)

#============================================================
# Input panel, enter details for wavefrom analysis

station = 'R7FA5' # Place Station name in this box
year = '2019'
month = '08'
day = '08'

#============================================================
# Client strictly for Raspberry Shake Stations

client = Client('http://fdsnws.raspberryshakedata.com')

# Set dt and dt1 to 30mins to an hours worth of data, 
# to allow the filtering to have a wider span of data

dt = UTCDateTime(str(year) + '-' + str(month) + '-' + str(day) + 'T16:30:00') 
dt1 = UTCDateTime(str(year) + '-' + str(month) + '-' + str(day) + 'T17:00:00')

#e1 and e2 are specific trimming tools, and are set for the known times of p and s wave arrivals.

e1 = UTCDateTime(str(year) + '-' + str(month) + '-' + str(day) + 'T16:52:05')
e2 = UTCDateTime(str(year) + '-' + str(month) + '-' + str(day) + 'T16:52:30')

#============================================================
# The try clause is used as the error codes on the client function are misleading
# If an error occurs, please insure the station name is correct.
# If it still throws and error, ensure the station is online during the time period you've selected, from IRIS
# I've found that even if some stations say they're online, sometimes there is no data to be streamed

try:
    st = client.get_waveforms('AM', station, '00', 'EHZ', dt, dt1) 
except:
    print('No information found from stations, please check that your entered values are correct, i.e station name and time')
    print('If entered data is correct, the station is not online for the selected time period')
else:
    
# Retrieving instrument response from station 
  
    inventory = client.get_stations(network = 'AM',
                                    station = station, 
                                    location = '00',
                                    channel = 'EHZ',
                                    level = 'response',
                                    starttime = dt,
                                    endtime = dt1)
           
    st_orig = st.copy()
    st_disp = st.copy()
 
# Applying pre filter
    
    pre_filt = (0.01, 0.02, 40, 45) 
#============================================================
# Poles and Zeros for WA simulation

WOODANDERSON = {'poles': [-5.49779 + 5.60886j, -5.49779 - 5.60886j],
        'zeros': [0 + 0j, 0 + 0j], 'gain': 1.0, 'sensitivity': 2080}

#============================================================
# Setting up Wood - Anderson profiles

st_WA = st.copy()

st_WA.remove_response(inventory = inventory, 
                        output = 'DISP', 
                        pre_filt = pre_filt)

tr_orig = st_orig[0]
tr_WA = st_WA[0]

# Simulating W-A instrument with known PaZ
tr_WA.simulate(paz_simulate = WOODANDERSON,
               simulate_sensitivity = False)

#============================================================
# Post filter and Trim

tr_orig.filter(type = 'highpass', freq = 1.25, corners = 4)
tr_WA.filter(type = 'highpass', freq = 1.25, corners = 4)

#Trimming to selected view period

tr_orig.trim(e1,e2)
tr_WA.trim(e1, e2)

# Calculating time x ticks from data.
# This is plot elapsed time from the beginning of the trim.
 
time = tr_orig.times()
#============================================================
# Settings for plots

f = plt.figure(figsize = (10,10))

ax = f.add_subplot(211)
ax2 = f.add_subplot(212)

ax.plot(time, tr_orig.data, 'k')
ax.set_ylabel(str(station) + ' [counts]')

ax2.plot(time, tr_WA.data, 'k')
ax2.set_ylabel('WA Displacment [m]')
ax2.set_xlabel('Time [s]')

plt.show()
#============================================================
#============================================================
