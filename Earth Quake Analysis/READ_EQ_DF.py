#! /usr/bin/env python3.7

import pandas as pd
from obspy.core.utcdatetime import UTCDateTime
from obspy.clients.fdsn import Client
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import seaborn as sns

para = {'axes.labelsize': 20, 'font.size': 20, 'legend.fontsize': 11,
'xtick.labelsize': 20,'ytick.labelsize': 20, 'figure.subplot.left': 0.12,
'figure.subplot.right': 0.98, 'figure.subplot.bottom': 0.11,
'figure.subplot.top': 0.97,'axes.linewidth': 1.0}
plt.rcParams.update(para)
#==================================================================================

# Following on from the Waveform analysis this code reads the generated CSV file
# and plots histos/scatter plots for the data

#==================================================================================

def read_EQ_df(filename):
    df  = pd.read_csv(filename,dtype={'Trace': object})

    df[['networkcode','Times','TrInfo']] = df['Trace'].str.split('|',expand=True)
    df[['StartTime','EndTime']] = df['Times'].str.split(' - ',expand=True)

    return df

if __name__=='__main__':
  
  filename = 'EQ_ANALYSED_DATAFRAME_LIB.csv'

  df = read_EQ_df(filename)
  # Creating Residual magnitude column
  df['ML BGS - ML Station Magnitude'] = df["BGS Magnitude"] - df["Station Magnitude"]

  pre_filt = (0.01, 0.02, 35, 45)
  client = Client('http://fdsnws.raspberryshakedata.com')

  

#==================================================================================
  # This allows the user to specifically targer one waveform to look at it's features
  # x is the index in the CSV file.
  x = 5
  
  stat = df.iloc[x,9]
  station = stat[3:8]
    
  
  StartTime = df.iloc[x,12]
  StartTime = UTCDateTime(StartTime)
  amp = df.iloc[x,6]
  
      #print(StartTime)
      
  EndTime = df.iloc[x,13]
  EndTime = UTCDateTime(EndTime)
      #print(EndTime)
  start = str(StartTime)
  end = str(EndTime)
      
  StartTime = UTCDateTime(StartTime - timedelta(seconds = 10))
  EndTime = UTCDateTime(EndTime + timedelta(seconds = 10))
  p = df.iloc[x, 7]
  s = df.iloc[x, 8]
  
  print('Amplitude :', amp)
  print('P-Wave :', p, ': S-Wave :', s)
  print('Trimmed Time : ', start[11:19], end[11:19])
  
      #print(StartTime, EndTime)
      
      
      
  res_magnitude = df.iloc[x, 14]
      #print(res_magnitude, ' ', station, ' ', StartTime, ' ', EndTime)
      
      
      
  st = client.get_waveforms('AM', station, '00', 'EHZ', StartTime, EndTime)
      
  inventory = client.get_stations(network = 'AM', station = station,
                                  location = '00',
                                  channel = 'EHZ',
                                  level = 'response',
                                  starttime = StartTime,
                                  endtime = EndTime)
            
  st.remove_response(inventory = inventory, output = 'DISP', pre_filt = pre_filt)
       
  tr = st[0]
      
  tr.filter(type = 'highpass', freq = 1.25, corners = 4)
      
  #tr.plot()
      
#==================================================================================
  
  res_mag = df['ML BGS - ML Station Magnitude']
  
  df2 = pd.DataFrame(df, columns =['Station Distance', 'ML BGS - ML Station Magnitude'])
  
  # This plots the scatter plot
  
  fig = sns.jointplot(x="Station Distance", 
                y="ML BGS - ML Station Magnitude", 
                data=df2, kind="scatter", 
                xlim =(0,400), ylim = (-1.5, 1.5))
  
  fig.ax_joint.set_xlabel('Station Distance [km]')
  
  plt.show(fig)

#==================================================================================
   
  sns.distplot(res_mag, 
             hist_kws={'edgecolor':'black', 'color' : 'darkslategrey' , 'alpha' : 0.8},
             kde= False,
             bins = (len(res_mag/0.2)))
  

  plt.ylabel('Number of Observations')
  plt.xlabel('ML BGS - ML Station Magnitude')
  plt.xlim(-1.5,1.5)   

  plt.plot()
  
  # A new csv file with more data.      
  df.to_csv('Picked_Data.csv', index=False)
  
      
      
      






