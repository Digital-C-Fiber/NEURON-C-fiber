import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plotRealAndSimulation(simulatedFilepath, saveFile, previousStim):
    data = pd.read_csv(simulatedFilepath,header=None, parse_dates=True)
    
    if previousStim:
        #drop first data points
        dropIndex=range(0,405)
        data=data.drop(dropIndex)
        data.reset_index(drop=True, inplace=True)
        for i, row in data.iterrows():
            data.loc[i,0] = data.loc[i,0]-360*1000

    for i, row in data.iterrows():
        data.loc[i,0] = data.loc[i,0]/1000

    
    dropIndex=[]
    for i, row in data.iterrows():
        if i+1< data.shape[0]:
            if data.loc[i+1,0]-data.loc[i,0]<3:
                dropIndex.append(i)

    #pd.set_option("display.max_rows", None, "display.max_columns", None)
    #print(data)
    data=data.drop(dropIndex)

    data_aps = np.genfromtxt('Protocols/Vivien/timestamps_aps.txt', delimiter="", names=["aps"])
    data_stim = np.genfromtxt('Protocols/Vivien/timestamps_stim.txt', delimiter="", names=["stim"])

    l = np.zeros((len(data_stim)-1,2))
    j=0
    i=0
    while i < len(data_stim):
        if j < len(data_aps):
            point = data_aps['aps'][j]-data_stim['stim'][i]
            if point>0.13:
                l[i][0]=data_aps['aps'][j]
                l[i][1]=(point*1000)
                j=j+1
                i=i+1
            else:
                j=j+1
        else:
            break

    data2 = pd.DataFrame(l)
    dropIndex=[]
    for i, row in data2.iterrows():
        if i+1< data2.shape[0]:
            if data2.loc[i+1,0]-data2.loc[i,0]<3:
                dropIndex.append(i)

    data2=data2.drop(dropIndex)
    #data2 = data2[data2.Latency != -1]
    #data2 = data2[data2.Latency>0.045]
    
    font = {'size': 15}
    plt.rc('font', **font)

    fig, ax1 = plt.subplots(1,1)
    data.plot.scatter(x=0, y=1, s=0.97,figsize=(15, 7.5), ax=ax1, c='red', secondary_y=True, label='Simulated Data')
    ax2 = ax1.twinx()
    data2.plot.scatter(x=0, y=1, s=0.97,figsize=(15, 7.5), ax=ax2, c='blue', secondary_y=True, label='Real Data')
    ax1.set_xlabel("time (s)")
    ax1.set_ylabel("Latency Simulated (ms)")
    ax2.set_ylabel("Latency Real (ms)")
    ax1.legend(loc=0)
    #ax2.legend(loc=2)
    ax1.axvline(l[342][0], color="red", linestyle="-")
    #ax1.set_ylim([6.8, 8.15])
    plt.show()
    fig = ax1.get_figure()
    fig.savefig(saveFile)