"""
pyVIRUS
library for modelling Corona pandenic data
2024-04-21
Georg Kaufmann
"""

import numpy as np
import matplotlib.pyplot as plt
import sys,os,csv,datetime
import libVIRUS


#================================#
def coronaReadDataGlobal(path='data/',timeEnd = '7/1/21'):
    """
    Corona virus
    Load global data files
    input:
      path - path to data sets (default)
    output:
      times    - array of days (string format m/d/yyyy)
      dataConfirmed  - number of confirmed infections (per country per day)
      dataRecovered  - number of recovered persons (per country per day)
      dataDeath      - number of death persons (per country per day)
      ndataUsed      - number of data points used (limited by choice of timeEnd)
    """
    inname = ['corona_timeseries_confirmed.csv','corona_timeseries_recovered.csv','corona_timeseries_deaths.csv']
    # open file and read with csv reader
    dataConfirmed  = list(csv.reader(open(path+inname[0], 'r', newline=''), delimiter=','))
    ndataConfirmed = len(dataConfirmed[0])
    dataRecovered  = list(csv.reader(open(path+inname[1], 'r', newline=''), delimiter=','))
    ndataRecovered = len(dataRecovered[0])
    dataDeath      = list(csv.reader(open(path+inname[2], 'r', newline=''), delimiter=','))
    ndataDeath     = len(dataDeath[0])
    # get datetime strings from first file (skip first 4 columns ...)
    times = dataConfirmed[0][4:ndataConfirmed]
    print('ndataConfirmed: ',ndataConfirmed)
    print('ndataRecovered: ',ndataRecovered)
    print('ndataDeath:     ',ndataDeath)
    print(times[-1])
    # define cut-off day for data set from timeEnd, get no of data points to use (ndataUsed)
    ndataUsed = 0
    read = True
    while(read):
        if (times[ndataUsed] != timeEnd):
            ndataUsed += 1
        else:
            read = False
    print('Last date used: ',times[ndataUsed],ndataUsed)
    return times,dataConfirmed,dataRecovered,dataDeath,ndataUsed
    

#================================#
def coronaPlotDataGlobal(times,dataConfirmed,dataRecovered,dataDeath,ndataUsed):
    """
    Corona virus
    Plot global data files
    input:
      times    - array of days (string format m/d/yyyy)
      dataConfirmed  - number of confirmed infections (per country per day)
      dataRecovered  - number of recovered persons (per country per day)
      dataDeath      - number of death persons (per country per day)
      ndataUsed      - number of data points used (limited by choice of timeEnd)
    output:
      (to file)
    """
    # define arrays
    time                = np.empty(0)
    cor_world_confirmed = np.zeros(ndataUsed-4,dtype='int')
    cor_world_death     = np.zeros(ndataUsed-4,dtype='int')
    cor_world_recovered = np.zeros(ndataUsed-4,dtype='int')
    # create datetime string
    for itimes in times[0:ndataUsed-4]:
        sep = itimes.split('/')
        d = int(sep[1])
        m = int(sep[0])
        y = int(sep[2])+2000
        time = np.append(time,datetime.date(y,m,d))

    # extract data for world
    for line in dataConfirmed:
        if (line[0] != 'Province/State'):
            for i in range(4,ndataUsed):
                cor_world_confirmed[i-4] = cor_world_confirmed[i-4] + int(float((line[i])))
    for line in dataDeath :
        if (line[0] != 'Province/State'):
            for i in range(4,ndataUsed):
                cor_world_death[i-4] = cor_world_death[i-4] + int(float((line[i])))
    for line in dataRecovered:
        if (line[0] != 'Province/State'):
            for i in range(4,ndataUsed):
                cor_world_recovered[i-4] = cor_world_recovered[i-4] + int(float((line[i])))

    print ("Infected:  %s %10i" % (time[ndataUsed-5],cor_world_confirmed[ndataUsed-5]))
    print ("Deaths:    %s %10i" % (time[ndataUsed-5],cor_world_death[ndataUsed-5]))
    print ("Recovered: %s %10i" % (time[ndataUsed-5],cor_world_recovered[ndataUsed-5]))
    # plot world data
    scale = 1e6
    plt.figure(figsize=(10.0,6.0))
    plt.title('Corona infections (World)')
    plt.xlabel('Time')
    plt.xticks(rotation=45)
    plt.ylabel('Mill. people')
    plt.plot(time,cor_world_confirmed/scale,linewidth='2',marker='o',markersize='2',color='blue',label='confirmed')
    plt.plot(time,cor_world_death/scale,linewidth='2',marker='o',markersize='2',color='red',label='deaths')
    plt.plot(time,cor_world_recovered/scale,linewidth='2',marker='o',markersize='2',color='green',label='recovered')
    plt.text(time[ndataUsed-5],cor_world_confirmed[ndataUsed-5]/scale,str(cor_world_confirmed[ndataUsed-5]),
             fontsize=12,color='blue',horizontalalignment='center',verticalalignment='bottom')
    plt.text(time[ndataUsed-5],cor_world_death[ndataUsed-5]/scale,str(cor_world_death[ndataUsed-5]),
             fontsize=12,color='red',horizontalalignment='center',verticalalignment='bottom')
    plt.text(time[ndataUsed-5],cor_world_recovered[ndataUsed-5]/scale,str(cor_world_recovered[ndataUsed-5]),
             fontsize=12,color='green',horizontalalignment='center',verticalalignment='bottom')
    plt.grid()
    plt.legend()
    return


#================================#
def coronaGetDataCountry(times,
                         dataConfirmed,dataRecovered,dataDeath,ndataUsed,
                         province='',country='Germany'):
    """
    Corona virus
    Extract country data
    input:
      times    - array of days (string format m/d/yyyy)
      dataConfirmed    - number of confirmed infections (per country per day)
      dataRecovered    - number of recovered persons (per country per day)
      dataDeath        - number of death persons (per country per day)
      ndataUsed        - number of data points used (limited by choice of timeEnd)
      province,country - province and country name (defaults)
    output:
      time             - dates as datetime string
      countryConfirmed - number of confirmed cases per day 
      countryRecovered - number of recovered cases per day
      countryDeath     - number of dead cases per day
      countryInvected  - number of infected cases per day (calculated)
    """
    # define arrays
    time          = np.empty(0)
    countryConfirmed = np.zeros(ndataUsed-4,dtype='int')
    countryDeath     = np.zeros(ndataUsed-4,dtype='int')
    countryRecovered = np.zeros(ndataUsed-4,dtype='int')
    countryInvected  = np.zeros(ndataUsed-4,dtype='int')

    # create datetime string
    for itimes in times[0:ndataUsed-4]:
        sep = itimes.split('/')
        d = int(sep[1])
        m = int(sep[0])
        y = int(sep[2])+2000
        time = np.append(time,datetime.date(y,m,d))
        
    # extract data for specific country
    for line in dataConfirmed:
        if (line[0] == province and line[1] == country):
            for i in range(4,ndataUsed):
                countryConfirmed[i-4] = countryConfirmed[i-4] + int(float((line[i])))
    for line in dataRecovered:
        if (line[0] == province and line[1] == country):
            for i in range(4,ndataUsed):
                countryRecovered[i-4] = countryRecovered[i-4] + int(float((line[i])))
    for line in dataDeath:
        if (line[0] == province and line[1] == country):
            for i in range(4,ndataUsed):
                countryDeath[i-4] = countryDeath[i-4] + int(float((line[i])))
    countryInvected = countryConfirmed - countryDeath - countryRecovered
    return time,countryConfirmed,countryRecovered,countryDeath,countryInvected


#================================#
def coronaPlotDataCountry(time,countryConfirmed,countryRecovered,countryDeath,ndataUsed,country='Germany'):
    """
    Corona virus
    Plot country data
    input:
      time             - dates as datetime string
      countryConfirmed - number of confirmed cases per day 
      countryRecovered - number of recovered cases per day
      countryDeath     - number of dead cases per day
      countryInvected  - number of infected cases per day (calculated)
      ndataUsed        - number of data points used (limited by choice of timeEnd)
    output:
      (to file)
    """
    scale = 1e6
    plt.figure(figsize=(12.0,6.0))
    plt.title('Corona infections ('+country+')')
    plt.xlabel('Time')
    plt.xticks(rotation=45)
    plt.ylabel('Mill. people')
    #plt.plot(time,cor_invected/scale,linewidth='2',color='red',alpha=0.5,label='invected')
    plt.plot(time,countryConfirmed/scale,linewidth='4',color='blue',label='confirmed')
    plt.plot(time,countryDeath/scale,linewidth='4',color='black',label='dead')
    plt.plot(time,countryRecovered/scale,linewidth='4',color='green',label='recovered')
    
    #plt.text(time[ndataUsed-5],cor_invected[ndataUsed-5]/scale,str(cor_invected[ndataUsed-5]),
    #         fontsize=12,color='red',horizontalalignment='left',verticalalignment='center')
    plt.text(time[ndataUsed-5],countryConfirmed[ndataUsed-5]/scale,str(countryConfirmed[ndataUsed-5]),
             fontsize=12,color='blue',horizontalalignment='left',verticalalignment='center')
    plt.text(time[ndataUsed-5],countryDeath[ndataUsed-5]/scale,str(countryDeath[ndataUsed-5]),
             fontsize=12,color='black',horizontalalignment='left',verticalalignment='center')
    plt.text(time[ndataUsed-5],countryRecovered[ndataUsed-5]/scale,str(countryRecovered[ndataUsed-5]),
             fontsize=12,color='green',horizontalalignment='left',verticalalignment='center')
    plt.grid()
    plt.legend()
    return


#================================#
def coronaPlotInfectedCountry(time,countyConfirmed,countryRecovered,countryDeath,countyInvected,ndataUsed,country='Germany'):
    """
    Corona virus
    Plot country infections
    input:
      time             - dates as datetime string
      countryConfirmed - number of confirmed cases per day 
      countryRecovered - number of recovered cases per day
      countryDeath     - number of dead cases per day
      countryInvected  - number of infected cases per day (calculated)
    output:
      (to file)
    """
    scale = 1e6
    plt.figure(figsize=(12.0,6.0))
    plt.title('Corona infections ('+country+')')
    plt.xlabel('Time')
    plt.xticks(rotation=45)
    plt.ylabel('Mill. people')
    plt.ylim([0,1.0])
    plt.fill_between(time,countyInvected/scale,0,edgecolor=(1,0,0,1.0),facecolor=(1,0,0,0.3),alpha=0.5,label='infected')
    plt.grid()
    plt.legend()
    return


#================================#
#================================#
#================================#
#================================#
#================================#