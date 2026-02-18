import pandas as pd
import numpy as np
import os

#filename can't be too long, full path can't be more than 255 characters
#therefore values are rounded!
def getFileSuffix(prot=1, scalingFactor=1, tempBranch=32, tempParent=37, 
                gPump=-0.0047891, gNav17=0.10664, gNav17Parent=0.10664, gNav18=0.24271, gNav18Parent=0.24271, gNav19=9.4779e-05, 
               gKs=0.0069733, gKf=0.012756, gH=0.0025377, gKdr=0.018002, gKna=0.00042, vRest=-55,
               sine=False, ampSine=0.1, particleNr=0):
    '''
    r=6
    fileSuffix=('_Prot'+str(prot)
                    +'_scale'+str(scalingFactor)
                    +'_gPump'+str(round(gPump,r))
                    +'_gNav17'+str(round(gNav17,r))
                    +'_gNav17P'+str(round(gNav17Parent,r))
                    +'_gNav18'+str(round(gNav18,r))
                    +'_gNav18P'+str(round(gNav18Parent,r))
                    +'_gNav19'+str(round(gNav19,r))
                    +'_gKs'+str(round(gKs,r))
                    +'_gKf'+str(round(gKf,r))
                    +'_gH'+str(round(gH,r))
                    +'_gKdr'+str(round(gKdr,r))
                    +'_gKna'+str(round(gKna,r))
                    +'_vRest'+str(vRest)
                    #+'_sine'+str(sine)
                    #+'_ampSine'+str(ampSine)
                    +'.csv')
    '''
    fileSuffix=('_Prot'+str(prot)
                    +'_particle'+str(particleNr)
                    +'.csv')
    return fileSuffix

#filetype can be "potential", "spikes" or "stim"
def getFilename(path="Results", filetype="potential", prot=1, scalingFactor=1, tempBranch=32, tempParent=37, 
                gPump=-0.0047891, gNav17=0.10664, gNav17Parent=0.10664, gNav18=0.24271, gNav18Parent=0.24271, gNav19=9.4779e-05, 
               gKs=0.0069733, gKf=0.012756, gH=0.0025377, gKdr=0.018002, gKna=0.00042, vRest=-55,
               sine=False, ampSine=0.1, nr=0):
    #old
    '''
    fileSuffix=('_Prot'+str(prot)+'_scalingFactor'+str(scalingFactor)
                +'_tempBranch'+str(tempBranch)+'_tempParent'+str(tempParent)
                +'_gPump'+str(gPump)+'_gNav17'+str(gNav17)+'_gNav18'+str(gNav18)+'_gNav19'+str(gNav19)
                +'_gKs'+str(gKs)+'_gKf'+str(gKf)+'_gH'+str(gH)+'_gKdr'+str(gKdr)+'_gKna'+str(gKna)+'_vRest'+str(vRest)+'.csv')
    
    fileSuffix=('_Prot'+str(prot)
                +'_gPump'+str(gPump)+'_gNav17'+str(gNav17)+'_gNav18'+str(gNav18)+'_gNav19'+str(gNav19)
                +'_gKs'+str(gKs)+'_gKf'+str(gKf)+'_gH'+str(gH)+'_gKdr'+str(gKdr)+'_gKna'+str(gKna)+'_vRest'+str(vRest)
                +'_sine'+str(sine)+'_ampSine'+str(ampSine)+'.csv')
    '''
    #print("GetFilename")
    if filetype=="stim":
        fileSuffix = '_Prot'+str(prot)+'.csv'
    else:
        fileSuffix = getFileSuffix(prot, scalingFactor, tempBranch, tempParent, 
                gPump, gNav17, gNav17Parent, gNav18, gNav18Parent, gNav19, 
               gKs, gKf, gH, gKdr, gKna, vRest,
               sine, ampSine, nr)
    
    #fileSuffix="particle"+str(nr)+".csv"
    
    filename = path+'/'+filetype+fileSuffix
    return filename

#filetype can be "potential", "spikes" or "stim"
def getData(path="Results", filetype="potential", prot=1, scalingFactor=1, tempBranch=32, tempParent=37, 
        gPump=-0.0047891, gNav17=0.10664, gNav17Parent=0.10664, gNav18=0.24271, gNav18Parent=0.24271, gNav19=9.4779e-05, 
        gKs=0.0069733, gKf=0.012756, gH=0.0025377, gKdr=0.018002, gKna=0.00042, vRest=-55,
        sine=False, ampSine=0.1, nr=0):
    #print("getData")
    filename=getFilename(path, filetype, prot, scalingFactor, tempBranch, tempParent, 
        gPump, gNav17, gNav17Parent, gNav18, gNav18Parent, gNav19, 
        gKs, gKf, gH, gKdr, gKna, vRest,
        sine, ampSine, nr)
    #print(filename)
    if os.path.exists(filename):
        data = pd.read_csv(filename)#, index_col=0)#, parse_dates=True, header=None)
        return data
    return None


def calculateLatency(data_aps, data_stim, norm=True):
    l = np.zeros((len(data_stim),2))
    j=0
    i=0
    k=0
    while i < len(data_stim):
        if j < len(data_aps):
            point = data_aps["Axon 3 1"][j]-data_stim["StimTime"][i]
            #print(point)
            if not np.isnan(point):
                if point>0 and point < 250:#normal AP
                    l[k][0]=data_stim["StimTime"][i]/1000
                    l[k][1]=point
                    j=j+1
                    i=i+1
                    k=k+1
                elif point <0:#AP belongs to next stimulation pulse
                    j=j+1
                elif point>250:#AP does not exist
                    l[k][0]=data_stim["StimTime"][i]/1000
                    l[k][1]=float('Nan')
                    i=i+1
                    k=k+1
            else:#AP does not exist
                l[k][0]=data_stim["StimTime"][i]/1000
                l[k][1]=float('Nan')
                i=i+1
                k=k+1     
        else:#no more APs
            l[k][0]=data_stim["StimTime"][i]/1000
            l[k][1]=float('Nan')
            i=i+1
            k=k+1
    #print(l)
    #normalize
    if norm:
        first=l[0][1]
        j=0
        for i in l:
            l[j][1]=(i[1]/first-1)*100
            j=j+1
    return l

#calculate latency for long model
def calculateLatencyLong(data_aps, data_stim, norm=True):
    l = np.zeros((len(data_stim),2))
    j=0
    i=0
    k=0
    while i < len(data_stim):
        if j < len(data_aps):
            point = data_aps["Axon 3 1"][j]-data_stim["StimTime"][i]
            #print(point)
            if not np.isnan(point):
                if point>0 and point < 300:#normal AP
                    l[k][0]=data_stim["StimTime"][i]/1000
                    l[k][1]=point
                    j=j+1
                    i=i+1
                    k=k+1
                elif point <0:#AP belongs to next stimulation pulse
                    j=j+1
                elif point>300:#AP does not exist
                    l[k][0]=data_stim["StimTime"][i]/1000
                    l[k][1]=float('Nan')
                    i=i+1
                    k=k+1
            else:#AP does not exist
                l[k][0]=data_stim["StimTime"][i]/1000
                l[k][1]=float('Nan')
                i=i+1
                k=k+1     
        else:#no more APs
            l[k][0]=data_stim["StimTime"][i]/1000
            l[k][1]=float('Nan')
            i=i+1
            k=k+1
    #print(l)
    #normalize
    if norm:
        first=l[0][1]
        j=0
        for i in l:
            l[j][1]=(i[1]/first-1)*100
            j=j+1
    return l

def calculateLatency2(data_aps, data_stim, norm=True):
    l = np.zeros((len(data_stim),2))
    j=0
    i=0
    k=0
    while i < len(data_stim):
        if j < len(data_aps):
            point = data_aps[j]-data_stim[i]
            #print(point)
            if point>0 and point < 400:
                l[k][0]=data_stim[i]/1000
                l[k][1]=point
                j=j+1
                i=i+1
                k=k+1
            elif point <0:
                j=j+1
            
            elif point>400:
                l[k][0]=data_stim[i]/1000
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

#for ELID and Serra
def getRealData(filename):
    
    data2 = pd.read_excel(filename, index_col=None, header=4)  

    data2 =data2.drop(1)
    data2 =data2.drop(labels="fiber id",axis=1)

    CM=[]
    CMi=[]
    VHT=[]
    for column in data2:
        #sort cm/cmi
        if data2.at[0,column] == "cm":
            CM.append(column)
        elif data2.at[0,column] == "cmi":
            CMi.append(column)
        elif data2.at[0,column] == "vht":
            VHT.append(column)

    data2= data2.drop(0)

    for column in data2:
        #normalize
        first=data2.at[2,column]
        for row in data2.index:
            data2.at[row,column]=(data2.at[row,column]/first-1)*100

    return data2, CM, CMi, VHT
    
    
def getRealDataODP(filename):
    data = pd.read_excel(filename, index_col=None, header=0)  

    CM=[]
    CMi=[]
    VHT=[]
    for index, row in data.iterrows():
        #sort cm/cmi
        if row['unit class'] == "cm":
            CM.append(index)
        elif row['unit class'] == "cmi":
            CMi.append(index)
        elif row['unit class'] == "vht":
            VHT.append(index)
    
    data=data.drop(['unit class', 'unit type'], axis=1)
    return data, CM, CMi, VHT
    