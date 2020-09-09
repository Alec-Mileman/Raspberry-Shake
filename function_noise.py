from obspy.core.utcdatetime import UTCDateTime
from obspy import read, read_inventory
from obspy.clients.fdsn import Client
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy as ob
from obspy.io.xseed import Parser
from datetime import datetime, timedelta

#======================================================================  

def time_selection(y, m, d):
    dt_s = [UTCDateTime(str(y) + '-' + str(m) + '-' + str(d) + 'T02:00:00'),
        UTCDateTime(str(y) + '-' + str(m) + '-' + str(d) + 'T08:00:00'),
        UTCDateTime(str(y) + '-' + str(m) + '-' + str(d) + 'T14:00:00'),
        UTCDateTime(str(y) + '-' + str(m) + '-' + str(d) + 'T20:00:00')]
    
    
    dt_e = [UTCDateTime(str(y) + '-' + str(m) + '-' + str(d) + 'T03:00:00'),
        UTCDateTime(str(y) + '-' + str(m) + '-' + str(d) + 'T09:00:00'),
        UTCDateTime(str(y) + '-' + str(m) + '-' + str(d) + 'T15:00:00'),
        UTCDateTime(str(y) + '-' + str(m) + '-' + str(d) + 'T21:00:00')]
    
    return dt_s, dt_e


def dir_read(day_number):
    file_name = {}
    file_name['E'] = '/home/alec/RSPSHK_DATA/2020/AM/RC54F/EHE.D/AM.RC54F.00.EHE.D.2020.'  + str(day_number)
    file_name['Z'] = '/home/alec/RSPSHK_DATA/2020/AM/RC54F/EHZ.D/AM.RC54F.00.EHZ.D.2020.'  + str(day_number)
    file_name['N'] = '/home/alec/RSPSHK_DATA/2020/AM/RC54F/EHN.D/AM.RC54F.00.EHN.D.2020.'  + str(day_number)
    
    return file_name 





def dir_read(day_number):
    
    file_name = '/home/alec/RSPBOOM_DATA/2020/AM/RDF08/HDF.D/AM.RDF08.00.HDF.D.2020.'  + str(day_number)
    
    return(file_name)

#======================================================================  

def dir_read_bkn(y,m,d):
        
    dir_bkn_mseed = []
    dir_bkn_mseed.append('/home/alec/BKNdata/BKN_BHZ_' + str(y) + '-' + str(m) + '-' + str(d) + 'T02:00:00.000000.mseed')
    dir_bkn_mseed.append('/home/alec/BKNdata/BKN_BHZ_' + str(y) + '-' + str(m) + '-' + str(d) + 'T08:00:00.000000.mseed')
    dir_bkn_mseed.append('/home/alec/BKNdata/BKN_BHZ_' + str(y) + '-' + str(m) + '-' + str(d) + 'T14:00:00.000000.mseed')
    dir_bkn_mseed.append('/home/alec/BKNdata/BKN_BHZ_' + str(y) + '-' + str(m) + '-' + str(d) + 'T20:00:00.000000.mseed')
    
    
    #dir_bkn_dlseed = []
    #dir_bkn_dlseed.append('/home/alec/BKNdata/BDSmetadata/BKN_BHZ_' + str(y) + '-' + str(m) + '-' + str(d) + 'T02:00:00.000000.dlseed')
    #dir_bkn_dlseed.append('/home/alec/BKNdata/BDSmetadata/BKN_BHZ_' + str(y) + '-' + str(m) + '-' + str(d) + 'T08:00:00.000000.dlseed')
    #dir_bkn_dlseed.append('/home/alec/BKNdata/BDSmetadata/BKN_BHZ_' + str(y) + '-' + str(m) + '-' + str(d) + 'T14:00:00.000000.dlseed')
    #dir_bkn_dlseed.append('/home/alec/BKNdata/BDSmetadata/BKN_BHZ_' + str(y) + '-' + str(m) + '-' + str(d) + 'T20:00:00.000000.dlseed')
    
    
    
    
    return dir_bkn_mseed#, dir_bkn_dlseed

#======================================================================  

def read_rm_resp(mseedpath,dlseedpath,resp='ACC'):
    '''
    mseedpath = miniseedpath
    dlseedpath = dataless seed path

    Returns obspy stream
    '''

    st = ob.read(mseedpath)

    # Used as a check for response removal. Ok.
    #st_orig = st.copy()
    #st_orig.detrend('linear')

    #---------- Remove instrument response -----------------------
    parser = Parser(dlseedpath)
    #print(parser.get_inventory())
    coords = parser.get_coordinates(st[0].stats.channel)
    paz = parser.get_paz(st[0].stats.channel)

    # BDS PAZ are already in DISP
    prefilt = (0.01,0.02,35,40)
    
    st.simulate(paz_remove=paz,zero_mean=True,pre_filt=prefilt)

    if resp=='ACC':
        st.differentiate()
        st.differentiate()
    elif resp=='VEL':
        st.differentiate()

    return st

#======================================================================  

def noise_prof_rspshk(day_number, year, month, day):
    
    day_num = day_number
    file_name = dir_read(day_num)
    pre_filt = (0.01, 0.02, 40.0, 45.0)
    
    #====================================================================
    
    
    dt_s, dt_e = time_selection(year, month, day)
    
    inv_RSPSHK_Z = read_inventory('/home/alec/Code/code_/out4.response.restored-EHZ-plus-decimation-new', format = 'RESP')
    inv_RSPSHK_E = read_inventory('/home/alec/Code/CODE/out4.response.restored-EHE-plus-decimation-new', format = 'RESP')
    inv_RSPSHK_N = read_inventory('/home/alec/Code/CODE/out4.response.restored-EHN-plus-decimation-new', format = 'RESP')
    
    xx_E = []
    xx_N = []
    xx_Z = []
    
    # Automate date selection!!
    
    for x in range(len(dt_s)):
        
            
            dt = dt_s[x]
            dt1 = dt_e[x]
    
            st_E = read(str(file_name['E']), startime = dt, endtime = dt1)
            st_Z = read(str(file_name['Z']), startime = dt, endtime = dt1)
            st_N = read(str(file_name['N']), startime = dt, endtime = dt1)
    
            
            st_E.remove_response(inventory = inv_RSPSHK_E, output = 'ACC', pre_filt = pre_filt)
            st_Z.remove_response(inventory = inv_RSPSHK_Z, output = 'ACC', pre_filt = pre_filt)
            st_N.remove_response(inventory = inv_RSPSHK_N, output = 'ACC', pre_filt = pre_filt)
    
            xx_E.append(st_E)
            xx_Z.append(st_Z)
            xx_N.append(st_N)
            
    return(xx_Z)

#======================================================================  

def noise_prof_rspshk_z(day_number, year, month, day):
    
    day_num = day_number
    file_name = dir_read(day_num)
    pre_filt = (0.01, 0.02, 40.0, 45.0)
    
    #====================================================================  
    dt_s, dt_e = time_selection(year, month, day)
    
    inv_RSPSHK_Z = read_inventory('/home/alec/Code/code_/out4.response.restored-EHZ-plus-decimation-new', format = 'RESP')

    xx_Z = []
    
    for x in range(len(dt_s)):
        
            
            dt = dt_s[x]
            dt1 = dt_e[x]
    
            st_Z = read(str(file_name['Z']), startime = dt, endtime = dt1)
            st_Z.remove_response(inventory = inv_RSPSHK_Z, output = 'ACC', pre_filt = pre_filt)
            print('stream', st_Z)
            xx_Z.append(st_Z)
            
            
    return(xx_Z)
#======================================================================  

def station_stream(y,m,d,pre_filt):
    
    client = Client('http://fdsnws.raspberryshakedata.com')

    station_num = ['R8118', 'R453F', 'R98E6', 'RFB08', 'R0353']
    
    dt_s, dt_e = time_selection(y,m,d)       
    
    at = []
    inv = []

    # Iterating the loop for the 4 selected time zones
    for x in range(len(dt_s)):
        
        dt = dt_s[x]
        dt1 = dt_e[x]
        
        # Response removal for all 4 stations
        for i in range(len(station_num)):
             
            st = client.get_waveforms('AM', station_num[i], '00', 'EHZ', dt, dt1)
            
            
            inventory = client.get_stations(network = 'AM',
                                            station = station_num[i],
                                            location = '00',
                                            channel = 'EHZ',
                                            level = 'response',
                                            starttime = dt,
                                            endtime = dt1)
            inv.append(inventory)
                   
            
            st.remove_response(inventory = inv[i],
                            output = 'ACC',
                            pre_filt = pre_filt)
            
            at.append(st)


    return(at)

def station_stream_RSPBOOM(y,m,d,pre_filt):
    
    client = Client('http://fdsnws.raspberryshakedata.com')

    station_num = ['R5754', 'R7BC1', 'RD120', 'RD5DB']
    
    dt_s, dt_e = time_selection(y,m,d)       
    
    at = []
    inv = []

    # Iterating the loop for the 4 selected time zones
    for x in range(len(dt_s)):
        
        dt = dt_s[x]
        dt1 = dt_e[x]
        
        # Response removal for all 4 stations
        for i in range(len(station_num)):
             
            st = client.get_waveforms('AM', station_num[i], '00', 'HDF', dt, dt1)
            
            
            inventory = client.get_stations(network = 'AM',
                                            station = station_num[i],
                                            location = '00',
                                            channel = 'HDF',
                                            level = 'response',
                                            starttime = dt,
                                            endtime = dt1)
            inv.append(inventory)
                   
            
            st.remove_response(inventory = inv[i],
                            output = 'ACC',
                            pre_filt = pre_filt)
            
            at.append(st)


    return(at)

#======================================================================  
def day_number(y,m,d):
    import datetime
    
    date = datetime.datetime.strptime(str(y) + '-' + str(m) + '-' + str(d), "%Y-%m-%d")
    
    x = date.timetuple()
    num = x.tm_yday
    
    return num
    
    

#pre_filt = (0.01, 0.02, 40.0, 45.0)
#abcd = station_stream('2020', '07', '21', pre_filt)
#abc = noise_prof_rspshk_z('203', '2020', '07', '21')
