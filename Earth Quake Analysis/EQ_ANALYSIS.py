import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn import Client
from geographiclib.geodesic import Geodesic
from datetime import datetime, timedelta
import math as m 

#=====================================================================
# This code reads from two text files, presenting EQ waveforms and allowing 
# the user to pick specific waves.
#
#=====================================================================

def add_EQ_dist(df_stations,df_eq,eq_index):
    EQlat = df_eq.iloc[eq_index]['LATITUDE']
    EQlon = df_eq.iloc[eq_index]['LONGITUDE']
    EQdate = df_eq.iloc[eq_index]['DATE']

    # Filter the stations, to provide only those online at time of EQ
    
    mask = (df_stations['StartTime'] < EQdate) & (df_stations['EndTime'] >EQdate)
    
    df_stations = df_stations.loc[mask]
    df_stations['Distance'] = np.zeros(len(df_stations))

    for row in df_stations.itertuples():
        azirge = Geodesic.WGS84.Inverse(EQlat,
                                        EQlon,
                                        row.Latitude,
                                        row.Longitude)
        dist = azirge['s12']/1000.
        df_stations.at[row.Index,'Distance'] = dist

        df_stations = df_stations.sort_values(by='Distance')

    return df_stations

def find_closest_stats(statfile,EQfile,eq_index):

    #------- Read in stations

    df_stations = pd.read_csv(statfile,delimiter='|',
                names=['Network', 'Station', 'Latitude', 'Longitude', 'Elevation',
                'Sitename', 'StartTime', 'EndTime'],
                skip_blank_lines=True,comment='#'
                )

    df_stations['StartTime'] = pd.to_datetime(df_stations['StartTime'])
    df_stations.EndTime = df_stations.EndTime.str.replace('2599','2050') # Pandas ran out of time
    df_stations['EndTime'] = pd.to_datetime(df_stations['EndTime'])



    #------ Read in the EQs

    df_eq = pd.read_csv(EQfile,delimiter=',',header=0)
    df_eq = df_eq.rename(columns=lambda x: x.strip()) # Remove whitespace from column names
    df_eq['DATE'] = pd.to_datetime(df_eq['DATE'])

    #------ Add distance for a particular EQ
    df_stations = add_EQ_dist(df_stations, df_eq, eq_index)

    return df_stations, df_eq.iloc[eq_index]

if __name__=='__main__':
    
        histo_mag = []
        dist_final = []
        eq_number = []
        trace_dataframe = []
        trace_stations = []
        mag_streams = []
        lat_dic = []
        lon_dic = []
        bgs_mag = []
        local_mag = []
        amplitude = []
        p_wave_b = []
        s_wave_b = []
        
        
        column_names = ['EQ Latitude', 'EQ Longitude',
                        'BGS Magnitude', 
                        'Station Distance',
                        'Station Magnitude', 'Trace', 
                        'Amplitude', 'P-Wave', 'S-Wave']
        
        df_eq_data = pd.DataFrame(columns = column_names)

   #===================================================================================================================================  
    
        for x in range(0,21,1): # Number of EQs in text file is 21. This loop iterates over each EQ
            
            stat_dist = []
            statfile = 'UK_EQ_STATIONS.txt'
            eqfile = 'UK_EQ.txt'
            df_stations, EQinfo = find_closest_stats(statfile, eqfile, eq_index = x) 
        
            print('-------')
            temp_eq = pd.DataFrame(EQinfo)  #Storing EQ info into a temporary dataframe
            print(EQinfo) # The Earthquake Parameters
            
    #===================================================================================================================================
        
            df_stations.reset_index(drop = True, inplace = True) # Resets index, ordering the closest station to epicenter as first index
            
            print('--o0o--')
            print('EQ ', x, ' out of 21')
    
    #===================================================================================================================================    
        # This section finds the amount of stations within 300km, adding the value to rnge. This is used later.
            
            temp_stations = pd.DataFrame(df_stations)
            temp_stations = temp_stations.applymap(str)
            rnge = 0
            
            for i in range(len(temp_stations)):
                distance = float(temp_stations.iloc[i, 8])
                
                if distance <= 300:
                    rnge += 1
    #=================================================================================================================================== 
            
            client = Client('http://fdsnws.raspberryshakedata.com') 
            station_mseed = []
    
            temp_eq = pd.DataFrame(EQinfo)
            temp_eq = temp_eq.applymap(str)
            
            t1 = temp_eq.loc['TIME'].values[0] # Grabs time of earthquake as a string
            t1 = t1[1:] # Removing Space in front of t1
            d1 = temp_eq.loc['DATE'].values[0] # Grabs date as a string
            d1 = d1[:10] # Removes time component of string
            
            dt = str(d1) + 'T' + str(t1) # Attatches date and time together in UTC datetime format
            dt = UTCDateTime(dt) # converts in UTCDateTime format for inventory read
            dt_eq = dt
            
            dt = dt - timedelta(minutes = 30) # dt is used to get an hours worth of data for filtering to work  
            dt1 = UTCDateTime(dt  + 3600) 
    
            s_wave = []
            p_wave = []
    
    #===================================================================================================================================      
    
            from UK_traveltimes import dist_to_arrtimes
            
            phaselist = ('p','s','P','S','Pg','Pn','Sg','Sn')
            pre_filt = (0.01, 0.02, 35, 45)
            streams = []
            WA_streams = []
            
            available = rnge
            stat_name = []
            
            arr_time = []
            end_time = []
            available = rnge
            
            WOODANDERSON = {'poles': [-5.49779 + 5.60886j, -5.49779 - 5.60886j],
                    'zeros': [0 + 0j, 0 + 0j], 'gain': 1.0, 'sensitivity': 2080}
            
    #===================================================================================================================================
            # This loops operates for all the stations within 300km, for each individual EQ.
            for j in range(rnge):
                    
                    progress = int((j/rnge)*100) # Progress updates for user
                    print(progress, '%')
                    
                    station = temp_stations.iloc[j, 1] # Getting station infor
                    
                    distance = float(temp_stations.iloc[j, 8])
                    
                    # Use of Try/Except clause as stations can be offline. This prevents script from breaking
                    try:
                        st = client.get_waveforms('AM', station, '00', 'EHZ', dt, dt1)
          
                    except:
                        
                        print('No data available at station:', station)
                        
                        available -=  1
                        
                        continue
                    
                    else: 
                        
                        depth = temp_eq.loc['DEPTH'].values[0]
                        depth = float(depth) 
                        stat_name.append(station)
                        stat_dist.append(distance)
                        
                        # Getting Inventory from IRIS to remove respone of station.
                        
                        inventory = client.get_stations(network = 'AM', station = station, 
                                                location = '00',
                                                channel = 'EHZ',
                                                level = 'response',
                                                starttime = dt,
                                                endtime = dt1)
                    
                       
                        # Calculating predicted arrival times depending on station dist and depth from EQ
                        
                        arrtimes, mintime, maxtime = dist_to_arrtimes(distance, depth, phaselist=phaselist)
                        
                        p_wave.append(float(mintime))
                        mintime = float(mintime) * 0.5 # Allowing a buffer of time for arrival
                        
                        s_wave.append(float(maxtime))
                        maxtime = float(maxtime) * 1.33
                        
                        e1 = UTCDateTime(dt_eq + timedelta(seconds = mintime))# Trimming of waveform for specific arrival time
                        
                        e2 = UTCDateTime(dt_eq + timedelta(seconds = maxtime)) 
                        
                        st_WA = st.copy()
                        
                        # Removing Response
                        st.remove_response(inventory = inventory, output = 'DISP', pre_filt = pre_filt)
                        
                        st_WA.remove_response(inventory = inventory, output = 'DISP', pre_filt = pre_filt)
                
                        tr = st[0]
                        tr_WA = st_WA[0]
                       
                        tr_WA.simulate(paz_simulate = WOODANDERSON, simulate_sensitivity = False)
                       
                        # Post Filter
                        tr.filter(type = 'highpass', freq = 1.25, corners = 4)
                        tr_WA.filter(type = 'highpass', freq = 1.25, corners = 4)
                      
                        tr.trim(e1, e2)
                        tr_WA.trim(e1, e2)

                        streams.append(tr)
                        WA_streams.append(tr_WA)
         
                                              
       
    #===================================================================================================================================
           
           # Notifying user of station availability 
            print(available,' out of', rnge, ' station streams available')
            perc = (available/rnge)*100
            print(int(perc), ' % available')

            
            # This section plots an image of the wavefrom to the user, allowing the to 'pick' the waves
            
            for k in range(len(WA_streams)):
                # Try/Except clause utilised as some streams were corrupt.
                try:
                    
                    WA_streams[k].plot()
                
                except:
                    print('Graph Error')
                    continue
                
                else:
                    
                    # Grabbing P and S Wave arrival times to display to user
                    eq_p = UTCDateTime(dt_eq + timedelta(seconds = p_wave[k]))
                    eq_p1 = str(eq_p)
                    eq_p1 = eq_p1[11:19]
                
                    eq_s = UTCDateTime(dt_eq + timedelta(seconds = s_wave[k]))
                    eq_s1 = str(eq_s)
                    eq_s1 = eq_s1[11:19]
                    
                    
                    # Asking user whether the waveform has clear/visible seismic signals
                
                    print('P-Wave:' + str(eq_p1) + '\n  S-Wave:' + str(eq_s1) + '\n')
                    
                    f = input('Is this graph suitable? [y/n]:')
                    
                    # Basic error catching implemented, incase the user inputs something other than [y/n]
                    
                    while f != 'y' and f != 'n':
                        
                        f = input('Please input correct value [y/n]')
                        
                        if f == 'y' or f =='n':
                            break
                
                if f == "y":
                    
                    # Creating arrays to store information about stations/waveforms
                    # to be used in analysis later/ csv files
                    
                    name = WA_streams[k].get_id()
                    
                    p_wave_b.append(eq_p1)
                    s_wave_b.append(eq_s1)
                    
                    trace_stations.append(name[3:8])
                    
                    
                    dist_final.append(stat_dist[k]) # Appends the station distance to a list
                    
                    
                    lat = temp_eq.loc['LATITUDE'].values[0]
                    lat_dic.append(lat)
                    
                    lon = temp_eq.loc['LONGITUDE'].values[0]
                    lon_dic.append(lon)
                    
                    bgs = temp_eq.loc['MAGNITUDE'].values[0]
                    bgs_mag.append(bgs)
              
                    eq_number.append(x)
                    
                    # This next section asks whether the user can see either P or S waves
                    
                    p = input('What wave do you want to trim for? [p/s]:')
                    
                    while p != 'p' and p != 's':
                        
                        p = input('Please enter correct variable [p/s]')
                        
                        if p =='y' or p =='n':
                            break
                    
                    if p == 'p':
                        
                        # Asking the user to trim the plot further for amplitude analysis
                        while True:
                            try:
                                a = float(input('Enter time reduction (seconds):'))
                            except ValueError:
                                print('Incorrect Entry')
                                continue
                            else:
                                break
                        
                        e1 = UTCDateTime(eq_p - timedelta(seconds = a))

                        while True:
                            
                            try:
                                b = float(input('Enter time extension (seconds):'))
                            
                            except ValueError:
                                print('Incorrect Entry')
                                continue
                            
                            else:
                                break
                        
                        e2 = UTCDateTime(eq_p + timedelta(seconds = b))
                        WA_streams[k].trim(e1, e2)
                        mag_streams.append(WA_streams[k])
                        
                        
                        
                    elif p == 's':
                        
                        while True:
                            
                            try:
                                a = float(input('Enter time reduction (seconds):'))
                            
                            except ValueError:
                                print('Incorrect Entry')
                                continue
                            
                            else:
                                break
                        
                        e1 = UTCDateTime(eq_s - timedelta(seconds = a))
                       
                        while True:
                            try:
                                b = float(input('Enter time extension (seconds):'))
                            except ValueError:
                                print('Incorrect Entry')
                                continue
                            else:
                                break
                            
                        e2 = UTCDateTime(eq_s + timedelta(seconds = b))
                        WA_streams[k].trim(e1, e2)
                        
                        mag_streams.append(WA_streams[k])
                        
                    
                    elif f == 'n':
                        continue
                    
                    
    #===================================================================================================================================        
        # Analysing the selected and trimmed waveforms from earlier.
        for i in range(len(mag_streams)):
            R = dist_final[i] 
                
            try:
                amplitude.append(mag_streams[i].max())
                
                
                A = abs(mag_streams[i].max() * 1e9) # Finding amplitude of waveform             
                
                local_mag.append(m.log10(A) + (0.95 * m.log10(R)) + (0.00183 * R) - 1.76) # Magnitude calculation
                
                histo_mag.append(float(bgs_mag[i]) - local_mag[i]) # Change mag ?
                
                
            except:
                    
                print('error')
                    
                pass
                
            else:
                        
                plt.hist(histo_mag, bins = 30)
                plt.show()
                
                data_eq = {'EQ Latitude': lat_dic[i], 'EQ Longitude': lon_dic[i],
                           'BGS Magnitude': bgs_mag[i],
                           'Station Distance' : dist_final[i] ,
                           'Station Magnitude' : local_mag[i],
                           'Trace' : mag_streams[i], 'Amplitude' : amplitude[i],
                           'P-Wave' : p_wave_b[i], 'S-Wave' : s_wave_b[i]}
                
            
            
            
            df_eq_data = df_eq_data.append(data_eq, ignore_index = True)            
            df_eq_data.to_csv('EQ_ANALYSED_DATAFRAME_.csv', index=False)            
            
    #===================================================================================================================================
                
            
            
    
        
        
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   