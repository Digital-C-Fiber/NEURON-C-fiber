from neuron import h
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
    
def plotLatency(data_aps, data_stim):
    l = np.zeros((len(data_stim),2))
    j=0
    i=0
    while i < len(data_stim):
        if j < len(data_aps):
            point = data_aps[j]-data_stim[i]
            if point>0:
                l[i][0]=data_stim[i]/1000
                l[i][1]=point
                j=j+1
                i=i+1
            else:
                j=j+1
        else:
            break
    if len(l)>0:
        plt.figure(figsize=(15,5))
        if len(data_stim)==1:
            plt.scatter(l[:,0],l[:,1])
        else:
            plt.scatter(l[:,0],l[:,1])
        plt.xlabel('time (s)')
        plt.ylabel('latency (ms)')
        plt.title('Latency')
        plt.show()
    return l[:,1]


def calculateLatency(data_aps, data_stim, norm=True):
    l = np.zeros((len(data_stim),2))
    j=0
    i=0
    k=0
    while i < len(data_stim):
        if j < len(data_aps):
            point = data_aps["Axon 3 1"][j]-data_stim["StimTime"][i]
            #print(point)
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
    #normalize
    if norm:
        first=l[0][1]
        j=0
        for i in l:
            l[j][1]=(i[1]/first-1)*100
            j=j+1
    return l

#plot data from previously saved file
def plotFromFile(prot, nseg):
    filename = 'data_Prot'+str(prot)+'_nseg'+str(nseg)+'.csv'
    data = np.genfromtxt(filename, delimiter=",", names=["t", "v_beg", "v_beg2", "v_mid", "v_end"])
    
    plt.figure(figsize=(15,5))
    plt.plot(data['t'],data['v_beg'])
    plt.xlabel('t (ms)')
    plt.ylabel('v (mV)')
    plt.title('axon(0)')
    plt.show()

    plt.figure(figsize=(15,5))
    plt.plot(data['t'],data['v_beg2'])
    plt.xlabel('t (ms)')
    plt.ylabel('v (mV)')
    plt.title('axon(0.1)')
    plt.show()

    plt.figure(figsize=(15,5))
    plt.plot(data['t'],data['v_mid'])
    plt.xlabel('t (ms)')
    plt.ylabel('v (mV)')
    plt.title('axon(0.5)')
    plt.show()

    plt.figure(figsize=(15,5))
    plt.plot(data['t'],data['v_end'])
    plt.xlabel('t (ms)')
    plt.ylabel('v (mV)')
    plt.title('axon(1)')
    plt.show()
    
def getFilename(prot=1, scalingFactor=1, tempBranch=32, tempParent=37, 
                gPump=-0.0047891, gNav17=0.10664, gNav18=0.24271, gNav19=9.4779e-05, 
               gKs=0.0069733, gKf=0.012756, gH=0.0025377, gKdr=0.018002, gKna=0.00042, vRest=-55,
               sine=False, ampSine=0.1):
    #old
    '''
    fileSuffix=('_Prot'+str(prot)+'_scalingFactor'+str(scalingFactor)
                +'_tempBranch'+str(tempBranch)+'_tempParent'+str(tempParent)
                +'_gPump'+str(gPump)+'_gNav17'+str(gNav17)+'_gNav18'+str(gNav18)+'_gNav19'+str(gNav19)
                +'_gKs'+str(gKs)+'_gKf'+str(gKf)+'_gH'+str(gH)+'_gKdr'+str(gKdr)+'_gKna'+str(gKna)+'_vRest'+str(vRest)+'.csv')
    '''
    '''
    fileSuffix=('_Prot'+str(prot)
                +'_gPump'+str(gPump)+'_gNav17'+str(gNav17)+'_gNav18'+str(gNav18)+'_gNav19'+str(gNav19)
                +'_gKs'+str(gKs)+'_gKf'+str(gKf)+'_gH'+str(gH)+'_gKdr'+str(gKdr)+'_gKna'+str(gKna)+'_vRest'+str(vRest)
                +'_sine'+str(sine)+'_ampSine'+str(ampSine)+'.csv')
    '''
    fileSuffix=('_Prot'+str(prot)
                +'_gPump'+str(round(gPump,4))+'_gNav17'+str(round(gNav17,4))+'_gNav18'+str(round(gNav18,4))+'_gNav19'+str(round(gNav19,4))
                +'_gKs'+str(round(gKs,4))+'_gKf'+str(round(gKf,4))+'_gH'+str(round(gH,4))+'_gKdr'+str(round(gKdr,4))+'_gKna'+str(round(gKna,4))+'_vRest'+str(round(vRest,4))
                +'.csv')
    return fileSuffix

#filetype can be "potential" or "spikes"
def getData(path="Results", filetype="potential", prot=1, scalingFactor=1, tempBranch=32, tempParent=37, 
        gPump=-0.0047891, gNav17=0.10664, gNav18=0.24271, gNav19=9.4779e-05, 
        gKs=0.0069733, gKf=0.012756, gH=0.0025377, gKdr=0.018002, gKna=0.00042, vRest=-55,
        sine=False, ampSine=0.1):
    #old
    '''
    fileSuffix=('_Prot'+str(prot)+'_scalingFactor'+str(scalingFactor)
                +'_tempBranch'+str(tempBranch)+'_tempParent'+str(tempParent)
                +'_gPump'+str(gPump)+'_gNav17'+str(gNav17)+'_gNav18'+str(gNav18)+'_gNav19'+str(gNav19)
                +'_gKs'+str(gKs)+'_gKf'+str(gKf)+'_gH'+str(gH)+'_gKdr'+str(gKdr)+'_gKna'+str(gKna)+'_vRest'+str(vRest)+'.csv')
    '''
    fileSuffix=('_Prot'+str(prot)
                +'_gPump'+str(gPump)+'_gNav17'+str(gNav17)+'_gNav18'+str(gNav18)+'_gNav19'+str(gNav19)
                +'_gKs'+str(gKs)+'_gKf'+str(gKf)+'_gH'+str(gH)+'_gKdr'+str(gKdr)+'_gKna'+str(gKna)+'_vRest'+str(vRest)
                +'_sine'+str(sine)+'_ampSine'+str(ampSine)+'.csv')
    filename = path+'/'+filetype+fileSuffix
    data = pd.read_csv(filename, index_col=None)#, parse_dates=True, header=None)
    return data
        