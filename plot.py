from neuron import h
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from ast import literal_eval
    
def getLatency(data_aps, data_stim):
    l = np.zeros((len(data_stim),2))
    j=0
    i=0
    k=0
    while i < len(data_stim):
        if j < len(data_aps):
            point = data_aps["Axon 3 1"][j]-data_stim["StimTime"][i]
            if point>0 and point < 400:
                l[k][0]=data_stim["StimTime"][i]/1000
                l[k][1]=point
                
                j=j+1
                i=i+1
                k=k+1
            elif point <0:
                j=j+1
            elif point>400:
                l[k][0]=data_stim["StimTime"][i]/1000
                l[k][1]=-1
                i=i+1
                k=k+1
                
        else:
            break
    return l

def plotLatency(data_aps, data_stim):
    l = getLatency(data_aps, data_stim)
       
    if len(l)>0:
        plt.figure(figsize=(15,5))
        plt.scatter(l[:,0],l[:,1])
        plt.xlabel('time (s)')
        plt.ylabel('latency (ms)')
        plt.title('Latency')
        plt.show()


def lower(ion):
    return ion-ion*0.5

def upper(ion):
    return ion+ion*0.5

#plot particle positions and scores
def plotParticles(dimension, num_particles, num_iterations, name):
    gPump=0.0047891 
    gNav17=0.10664 
    gNav18=0.24271 
    gNav19=9.4779e-05 
    gKs=0.0069733 
    gKf=0.012756 
    gH=0.0025377 
    gKdr=0.018002 
    gKna=0.00042
    
    channel_names=["Pump", "Nav17", "Nav18", "Nav19", "Ks", "Kf", "h", "Kdr", "Kna"]
    ions=[gPump, gNav17, gNav18, gNav19, gKs, gKf, gH, gKdr, gKna]
    
    '''
    colors=['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    '''
    
    #get particles 
    for j in range(dimension):
        for i in range(num_particles):
            filename='Results/'+str(name)+'particle'+str(i)+'.csv'
            data = pd.read_csv(filename)

            data['Position'] = data['Position'].apply(literal_eval)#make string to list
            
            #weights = np.arange(1, len(data["Position"])+1)
            #plt.scatter(data["Score"], data["Position"].apply(lambda x: x[j]), marker='.', c=weights, cmap=colors[i])
            
            plt.plot(data["Score"], data["Position"].apply(lambda x: x[j]), marker='.')
        
        plt.ylim(lower(ions[j]), upper(ions[j]))
        plt.xlim(0,2)
        plt.title(channel_names[j])
        plt.xlabel("Score")
        plt.ylabel("conduction")
        plt.show()


def plotRecoveryCycle(data_aps, data_stim):
    l = getLatency(data_aps, data_stim)
    recoveryCycle=[]
    extras=[2000, 1750, 1500, 1250, 1000, 750, 500, 250, 150, 100, 75, 50, 40, 30, 20, 10]#distance to next regular puls
    extrasFinal=[]
    
    numberOfExtra=len(extras)#total number of extra pulses
    numRegPulses=6#number of regular pulses in between extra pulses
    for i in range(numberOfExtra):
        l1=l[20+i*numRegPulses][1]#latency of first extra pulse
        l2=l[21+i*numRegPulses][1]#latency of next regular pulse
        if l1!=-1 and l2!=-1:
            recoveryCycle.append(l2-l1)
            extrasFinal.append(extras[i])
    print(recoveryCycle)
    print(extrasFinal)
    
    plt.plot(extrasFinal, recoveryCycle, marker=".")
    plt.xlabel("interspike interval")
    plt.ylabel("slowing/speeding")
    


                