from neuron import h
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from ast import literal_eval
import dataProcessing
import seaborn as sns
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import matplotlib
    

def plotLatency(data_aps, data_stim, norm=False, title="Latency"):
    l = dataProcessing.calculateLatency(data_aps, data_stim, norm)
       
    if len(l)>0:
        plt.figure(figsize=(15,5))
        plt.scatter(l[:,0],l[:,1])
        plt.xlabel('time (s)')
        if norm:
            plt.ylabel('latency (%)')
        else:
            plt.ylabel('latency (ms)')
        plt.title(title)
        plt.show()


#get lower bound for channel
def lower(ion):
    vRest=-55
    if ion==vRest:
        return vRest-5
    return 0#because negative values don't make sense
    '''
    if ion==gKdr or ion==gNav18:
        return ion-2*ion
    else:
        return ion-4*ion
    '''
    #return ion-(0.5*ion)
    #return -4.5#beale
    #return -2*math.pi#mishra bird
    #return -10#holder
    #return -512#eggholder
    #return -100#schaffer, easom


#get upper bound for channel
def upper(ion):
    gPump=0.0047891 
    gNav17=0.10664 
    gNav18=0.24271 
    gNav19=9.4779e-05 
    gKs=0.0069733 
    gKf=0.012756 
    gH=0.0025377 
    gKdr=0.018002 
    gKna=0.00042
    vRest=-55
    
    if ion==gKdr or ion==gNav18:
        return 2*ion
    elif ion==gH:
        return 6*gH
    elif ion==gPump:
        return 5*ion
    elif ion==vRest:
        return vRest+10
    else:
        return 4*ion
    
    #return ion+(0.5*ion)
    #return 4.5#beale
    #return 2*math.pi#mishra bird
    #return 10#holder
    #return 512#eggholder
    #return 100#schaffer, easom

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
            filename='Results/'+str(name)+'/'+str(name)+'particle'+str(i)+'.csv'
            data = pd.read_csv(filename)

            data['Position'] = data['Position'].apply(literal_eval)#make string to list
            
            #weights = np.arange(1, len(data["Position"])+1)
            #plt.scatter(data["Score"], data["Position"].apply(lambda x: x[j]), marker='.', c=weights, cmap=colors[i])
            
            #plt.plot(data["Score"], data["Position"].apply(lambda x: x[j]), marker='.')
            plt.plot(data["Position"].apply(lambda x: x[j]), data["Score"], marker='.')
        
        plt.xlim(lower(ions[j]), upper(ions[j]))
        #plt.xlim(0,2)
        plt.title(channel_names[j])
        plt.ylabel("Score")
        plt.xlabel("conduction")
        plt.show()

    
def plotParticles2D(dimension, num_particles, num_iterations, name, showInf, scoreCutoff=float('inf'), path="", infNumber=100):
    gPump=0.0047891 
    gNav17=0.10664 
    gNav17Parent=0.10664 
    gNav18=0.24271
    gNav18Parent=0.24271 
    gNav19=9.4779e-05 
    gKs=0.0069733 
    gKf=0.012756 
    gH=0.0025377 
    gKdr=0.018002 
    gKna=0.00042
    vRest=-55
    
    channel_names=["Pump", "Nav17", "Nav17Parent", "Nav18", "Nav18Parent", "Nav19", "Ks", "Kf", "h", "Kdr", "Kna", "RMP"]
    ions=[gPump, gNav17, gNav17Parent, gNav18, gNav18Parent, gNav19, gKs, gKf, gH, gKdr, gKna, vRest]
    
    filename='Results/'+str(path)+'/'+str(name)+'/'+str(name)+'_bestParticles.csv'
    dataBest = pd.read_csv(filename)
    dataBest['Position'] = dataBest['Position'].apply(literal_eval)#make string to list
    
    #get maxium value
    maxScore=0
    '''
    for i in range(num_particles):
        if len(allData[i]["Score"]) != 0:
            #convert inf to negative inf, then find max
            maxi = max(np.where(np.isinf(allData[i]["Score"]),-np.Inf,allData[i]["Score"]))
            if maxScore<maxi:
                maxScore=maxi
    '''            

    maxInf=0
    allData=[]
    for i in range(num_particles):
        filename='Results/'+str(path)+'/'+str(name)+'/'+str(name)+'particle'+str(i)+'.csv'
        data = pd.read_csv(filename)
        data['Position'] = data['Position'].apply(literal_eval)#make string to list
        data=data[0:num_iterations]
        
        if scoreCutoff!=float('inf'):
            data = data[data['Score'] < scoreCutoff]#only data with score < scoreCutoff
            data=data.reset_index()
        #get maxium value
        if len(data["Score"]) != 0:
            #convert inf to negative inf, then find max
            maxi = max(np.where(np.isinf(data["Score"]),-np.Inf,data["Score"]))
            if maxScore<maxi:
                maxScore=maxi
        maxInf=maxScore+infNumber
        if showInf:
            data.replace([np.inf, -np.inf], maxInf, inplace=True)
        else:
            data=data.replace([np.inf, -np.inf], np.nan).dropna()
        allData.append(data)
     
    #print(allData)
    figure, axis = plt.subplots(dimension, dimension, figsize=(25, 25))
    #get particles 

    # getting the original colormap using cm.get_cmap() function
    orig_map=plt.cm.get_cmap('YlGnBu')
    # reversing the original colormap using reversed() function
    reversed_map = orig_map.reversed()
    
    #get minimum value
    length=len(dataBest)
    minScore=dataBest["Score"][length-1]
    
    norm = plt.Normalize(minScore, maxScore)
    
    for p1 in range(dimension):
        for p2 in range(dimension):
            dataHist=[]
            for i in range(num_particles):
                if p2<p1:
                    im=axis[p1,p2].scatter(allData[i]["Position"].apply(lambda x: x[p2]), allData[i]["Position"].apply(lambda x: x[p1]), c=reversed_map(norm(allData[i]["Score"])), cmap=reversed_map, marker=".", s=3)
                    axis[p1,p2].plot(dataBest["Position"][num_iterations-1][p2], dataBest["Position"][num_iterations-1][p1], c="Red", zorder=100, markersize=10, marker="*")
                elif p1==p2:
                    dataHist.append(allData[i]["Position"].apply(lambda x: x[p1]))
            '''
            if p1==p2:
                axis[p1,p2].hist(dataHist)
                if ions[p2]==gKdr or ions[p2]==gNav18:
                    axis[p1,p2].set_xlim(0,ions[p2]*2)
                elif ions[p2]==vRest:
                    axis[p1,p2].set_xlim(vRest-5,vRest+5)
                else:
                    axis[p1,p2].set_xlim(0,ions[p2]*4)
            elif p2<p1:
                if ions[p1]==gKdr or ions[p1]==gNav18:
                    axis[p1,p2].set_ylim(0,ions[p1]*2)
                elif ions[p1]==vRest:
                    axis[p1,p2].set_ylim(vRest-5,vRest+5)
                else:
                    axis[p1,p2].set_ylim(0,ions[p1]*4)

                if ions[p2]==gKdr or ions[p2]==gNav18:
                    axis[p1,p2].set_xlim(0,ions[p2]*2)
                elif ions[p2]==vRest:
                    axis[p1,p2].set_xlim(vRest-5,vRest+5)
                else:
                    axis[p1,p2].set_xlim(0,ions[p2]*4)
            '''
            if p1==p2:
                axis[p1,p2].hist(dataHist)
                axis[p1,p2].set_xlim(lower(ions[p2]),upper(ions[p2]))
            elif p2<p1:
                axis[p1,p2].set_xlim(lower(ions[p2]),upper(ions[p2]))
                axis[p1,p2].set_ylim(lower(ions[p1]),upper(ions[p1]))
            
            axis[p1,0].set_ylabel(channel_names[p1])
            axis[dimension-1,p2].set_xlabel(channel_names[p2])
    
    
    if maxScore>scoreCutoff:
        maxScore=scoreCutoff
    print("Max score: ", maxScore)  
    
    print("Min score: ", minScore)
    if showInf:
        normalizer=Normalize(minScore, maxInf)
    else:
        normalizer=Normalize(minScore,maxScore)
    im=cm.ScalarMappable(norm=normalizer, cmap=reversed_map)
    figure.colorbar(im, ax=axis.ravel().tolist())
    
def plotOneParticle2D(dimension, num_particles, num_iterations, name, showInf, scoreCutoff=float('inf'), path="", maxInf=200):
    gPump=0.0047891 
    gNav17=0.10664 
    gNav17Parent=0.10664 
    gNav18=0.24271
    gNav18Parent=0.24271 
    gNav19=9.4779e-05 
    gKs=0.0069733 
    gKf=0.012756 
    gH=0.0025377 
    gKdr=0.018002 
    gKna=0.00042
    
    channel_names=["Nav1.7", "Nav1.8"]
    ions=[gNav17,gNav18]
    
    filename='Results/'+str(path)+'/'+str(name)+'/'+str(name)+'_bestParticles.csv'
    dataBest = pd.read_csv(filename)
    dataBest['Position'] = dataBest['Position'].apply(literal_eval)#make string to list
    
    #maxInf=200
    allData=[]
    for i in range(num_particles):
        filename='Results/'+str(path)+'/'+str(name)+'/'+str(name)+'particle'+str(i)+'.csv'
        data = pd.read_csv(filename)
        data['Position'] = data['Position'].apply(literal_eval)#make string to list
        data=data[0:num_iterations]
        
        if scoreCutoff!=float('inf'):
            data = data[data['Score'] < scoreCutoff]#only data with score < scoreCutoff
            data=data.reset_index()
        
        if showInf:
            data.replace([np.inf, -np.inf], maxInf, inplace=True)
        else:
            data=data.replace([np.inf, -np.inf], np.nan).dropna()
        allData.append(data)
     
    font = {'family' : 'arial',
        'size'   : 20}

    matplotlib.rc('font', **font)
    
    #print(allData)
    figure, axis = plt.subplots(1, 1)
    #get particles 
    
    
    # getting the original colormap using cm.get_cmap() function
    orig_map=plt.cm.get_cmap('YlGnBu')
    # reversing the original colormap using reversed() function
    reversed_map = orig_map.reversed()
        
    #get maxium value
    maxScore=0
    for i in range(num_particles):
        if len(allData[i]["Score"]) != 0:
            #convert inf to negative inf, then find max
            maxi = max(np.where(np.isinf(allData[i]["Score"]),-np.Inf,allData[i]["Score"]))
            if maxScore<maxi:
                maxScore=maxi
                
    #get minimum value
    length=len(dataBest)
    minScore=dataBest["Score"][length-1]
    
    norm = plt.Normalize(minScore, maxScore)
    

    for i in range(num_particles):
        im=axis.scatter(allData[i]["Position"].apply(lambda x: x[3]), allData[i]["Position"].apply(lambda x: x[1]), c=reversed_map(norm(allData[i]["Score"])), cmap=reversed_map, marker=".")

        axis.plot(dataBest["Position"][num_iterations-1][3], dataBest["Position"][num_iterations-1][1], c="Red", zorder=100, markersize=10, marker="*")

    axis.set_ylabel(channel_names[0])
    axis.set_xlabel(channel_names[1])
    
    
    if maxScore>scoreCutoff:
        maxScore=scoreCutoff
    print("Max score: ", maxScore)  
    
    print("Min score: ", minScore)
    if showInf:
        normalizer=Normalize(minScore,maxInf)
    else:
        normalizer=Normalize(minScore,maxScore)
    im=cm.ScalarMappable(norm=normalizer, cmap=reversed_map)
    figure.colorbar(im, ax=axis)
    plt.savefig("Figures/Optimization-Particle", bbox_inches='tight', dpi=300)
    
    #figure.subplots_adjust(right=0.8)
    #cbar_ax = figure.add_axes([0.85, 0.15, 0.05, 0.7])
    #figure.colorbar(im, cax=cbar_ax)
#plots 2d marginals of particles as a heatmap, and 1d marginals as histogram 
#dimension: number of dimensions of particle
#num_particles: number of particles
#num_iterations: number of iterations
#name: filename 
#replaceInf: replace infinite values with number
#scoreCutoff: only plot particles with scores below cutoff
#nBins: number of bins shown in heatmap
def plotParticles2DMarginal(dimension, num_particles, num_iterations, name, showInf=False, scoreCutoff=float('inf'), maxInf=200, nBins=10):
    gPump=0.0047891 
    gNav17Parent=0.10664 
    gNav17=0.10664 
    gNav18=0.24271
    gNav18Parent=0.24271 
    gNav19=9.4779e-05 
    gKs=0.0069733 
    gKf=0.012756 
    gH=0.0025377 
    gKdr=0.018002 
    gKna=0.00042
    vRest=-55
    
    channel_names=["Pump", "Nav17", "Nav17Parent", "Nav18", "Nav18Parent", "Nav19", "Ks", "Kf", "h", "Kdr", "Kna", "RMP"]
    ions=[gPump, gNav17, gNav17Parent, gNav18, gNav18Parent, gNav19, gKs, gKf, gH, gKdr, gKna, vRest]
    
    #get best particles
    filename='Results/'+str(name)+'/'+str(name)+'_bestParticles.csv'
    dataBest = pd.read_csv(filename)
    dataBest['Position'] = dataBest['Position'].apply(literal_eval)#make string to list
    
    #get all particles 
    allData=[]
    for i in range(num_particles):
        filename='Results/'+str(name)+'/'+str(name)+'particle'+str(i)+'.csv'
        data = pd.read_csv(filename)
        data['Position'] = data['Position'].apply(literal_eval)#make string to list
        data=data[0:num_iterations]
        
        if scoreCutoff!=float('inf'):
            data = data[data['Score'] < scoreCutoff]#only data with score < scoreCutoff
            data = data.reset_index()
                
        if showInf:
            data.replace([np.inf, -np.inf], maxInf, inplace=True)
        else:
            data=data.replace([np.inf, -np.inf], np.nan).dropna()
        allData.append(data)
    
    figure, axis = plt.subplots(dimension, dimension, figsize=(25, 25))            
    # getting the original colormap using cm.get_cmap() function
    orig_map=plt.cm.get_cmap('YlGnBu')
    # reversing the original colormap using reversed() function
    reversed_map = orig_map.reversed()
    
    #get maxium value
    maxScore=0
    for i in range(num_particles):
        if len(allData[i]["Score"]) != 0:
            #convert inf to negative inf, then find max
            maxi = max(np.where(np.isinf(allData[i]["Score"]),-np.Inf,allData[i]["Score"]))
            if maxScore<maxi:
                maxScore=maxi
                
    #get minimum value
    length=len(dataBest)
    minScore=dataBest["Score"][length-1]
    
    norm = plt.Normalize(minScore, maxScore)
    
    for p1 in range(dimension):
        for p2 in range(dimension):
            dataHist=[]
            
            if p2<p1:
                gridFinal = np.zeros_like(allData[0]["Score"], shape=(nBins+1, nBins+1))
                countFinal = np.zeros_like(allData[0]["Score"], shape=(nBins+1, nBins+1))

                #set step size
                if ions[p1]==gKdr or ions[p1]==gNav18:
                    xstep=ions[p1]*2/nBins
                else:
                    xstep=ions[p1]*4/nBins
                #print(xstep)

                if ions[p2]==gKdr or ions[p2]==gNav18:
                    ystep=ions[p2]*2/nBins
                else:
                    ystep=ions[p2]*4/nBins
                #print(ystep)
 
            for i in range(num_particles):
                if p2<p1:
                    x=allData[i]["Position"].apply(lambda x: x[p1])#Nav17
                    y=allData[i]["Position"].apply(lambda x: x[p2])#pump
                    #print(x)
                    #print(y)
                    
                    #calculate bins
                    xbin = ((x) // xstep).astype(int)
                    ybin = ((y) // ystep).astype(int) 
                    #print(xbin)
                    #print(ybin)
                    
                    grid = np.zeros_like(allData[i]["Score"], shape=(nBins+1, nBins+1))
                    count = np.zeros_like(allData[i]["Score"], shape=(nBins+1, nBins+1))
                    np.add.at(grid, (xbin, ybin), allData[i]["Score"])# add scores
                    #print(grid)
                    np.add.at(count, (xbin, ybin), 1)#number of values in each field
                    #print("Count")
                    #print(count)
                    
                    #print("Grid")
                    #print(grid)
                    gridFinal+=grid
                    countFinal+=count

                elif p1==p2:
                    dataHist.append(allData[i]["Position"].apply(lambda x: x[p1]))
                 
            if p2<p1:              
                gridFinal=gridFinal/countFinal#get mean score
                #print("GridFinal")
                #print(gridFinal)
                
                #print("CountFinal")
                #print(countFinal)
                
                #find max score
                #maxi = np.where(np.isinf(gridFinal),-np.Inf, gridFinal).argmax()
                #test=np.where(np.isinf(gridFinal),-np.Inf, gridFinal)
                #maxi = np.nanmax(test)
                #maxi = max(map(max, gridFinal))
                #if maxScore<maxi:
                #    maxScore=maxi

                axis[p1,p2].imshow(gridFinal, cmap=reversed_map, origin='lower')
                
                #set labels
                if ions[p1]==gKdr or ions[p1]==gNav18:
                    xstepL=ions[p1]*2/5
                else:
                    xstepL=ions[p1]*4/5
                #print(xstep)

                if ions[p2]==gKdr or ions[p2]==gNav18:
                    ystepL=ions[p2]*2/3
                else:
                    ystepL=ions[p2]*4/3
                
                labelX = np.array([round(xstepL*i,4) for i in range(1,6)])
                labelY = np.array([round(ystepL*i,3) for i in range(1,4)])

                axis[p1,p2].set_xticks(np.arange(0,nBins, nBins/3))
                axis[p1,p2].set_yticks(np.arange(0,nBins, nBins/5))
                axis[p1,p2].set_xticklabels(labels=labelY)
                axis[p1,p2].set_yticklabels(labels='')
                #axis[0,p2].set_xticklabels(labels=labelY)
                axis[p1,0].set_yticklabels(labels=labelX)

                #mark best position
                y=dataBest["Position"][num_iterations-1][p2]
                x=dataBest["Position"][num_iterations-1][p1]

                xbin = ((x) // xstep)
                ybin = ((y) // ystep)

                #axis[p1,p2].plot(dataBest["Position"][num_iterations-1][p2], dataBest["Position"][num_iterations-1][p1], c="Red", zorder=100, markersize=10, marker="*")
                axis[p1,p2].plot(ybin, xbin, c="Red", zorder=100, markersize=10, marker="*")
                
            if p1==p2:
                axis[p1,p2].hist(dataHist)
                    
            axis[p1,0].set_ylabel(channel_names[p1])
            axis[dimension-1,p2].set_xlabel(channel_names[p2])
    
    if maxScore>scoreCutoff:
        maxScore=scoreCutoff
    print("Max score: ", maxScore)  
    print("Min score: ", minScore)
    if showInf:
        normalizer=Normalize(minScore,maxInf)
    else:
        normalizer=Normalize(minScore,maxScore)
    im=cm.ScalarMappable(norm=normalizer, cmap=reversed_map)
    figure.colorbar(im, ax=axis.ravel().tolist())
    
    #figure.subplots_adjust(right=0.8)
    #cbar_ax = figure.add_axes([0.85, 0.15, 0.05, 0.7])
    #figure.colorbar(im, cax=cbar_ax)    
    
def plotParticles2Dindividual(dimension, num_particles, num_iterations, name):
    gPump=0.0047891 
    gNav17Parent=0.10664 
    gNav17=0.10664 
    gNav18=0.24271
    gNav18Parent=0.24271 
    gNav19=9.4779e-05 
    gKs=0.0069733 
    gKf=0.012756 
    gH=0.0025377 
    gKdr=0.018002 
    gKna=0.00042
    vRest=-55
    
    channel_names=["Pump", "Nav17", "Nav17Parent", "Nav18", "Nav18Parent", "Nav19", "Ks", "Kf", "h", "Kdr", "Kna", "RMP"]
    ions=[gPump, gNav17, gNav17Parent, gNav18, gNav18Parent, gNav19, gKs, gKf, gH, gKdr, gKna, vRest]
    
    '''
    colors=['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                      'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                      'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    '''
 
    #get particles 
    for i in range(num_particles):
        figure, axis = plt.subplots(dimension, dimension, figsize=(25, 25))
        for p1 in range(dimension):
            for p2 in range(dimension):
                dataHist=[]
                #i=1

                filename='Results/'+str(name)+'/'+str(name)+'particle'+str(i)+'.csv'
                data = pd.read_csv(filename)
                data['Position'] = data['Position'].apply(literal_eval)#make string to list
                data=data[0:num_iterations]

                data.replace([np.inf, -np.inf], 10, inplace=True)
                #print(data["Score"])
                if p2<p1:
                    #im=axis[p1,p2].scatter(data["Position"].apply(lambda x: x[p2]), data["Position"].apply(lambda x: x[p1]), marker='.', c=data["Score"], cmap="seismic")
                    im=axis[p1,p2].scatter(data["Position"].apply(lambda x: x[p2]), data["Position"].apply(lambda x: x[p1]), c=range(num_iterations), cmap="Blues")
                elif p1==p2:
                    #w=upper(ions[p1])/len(data["Score"])
                    #axis[p1,p2].bar(data["Position"].apply(lambda x: x[p1]), data["Score"], width=w, color='c')
                    #b=np.arange(upper(ions[p1]), step=upper(ions[p1])/10)
                    #print(b)
                    dataHist.append(data["Position"].apply(lambda x: x[p1]))
                

                if p1==p2:
                    #print(dataHist)
                    axis[p1,p2].hist(dataHist)


                #axis[p1,p2].set_xlim(lower(ions[p2]), upper(ions[p2]))
                #axis[p1,p2].set_ylim(lower(ions[p1]), upper(ions[p1]))
                #plt.xlim(0,2)
                axis[p1,0].set_ylabel(channel_names[p1])
                axis[dimension-1,p2].set_xlabel(channel_names[p2])
        figure.subplots_adjust(right=0.8)
        cbar_ax = figure.add_axes([0.85, 0.15, 0.05, 0.7])
        figure.colorbar(im, cax=cbar_ax)
        
        
def plotRecoveryCycle(data_aps, data_stim, norm):
    l = dataProcessing.calculateLatency(data_aps, data_stim, norm)
    recoveryCycle=[]
    #extras=[2000, 1750, 1500, 1250, 1000, 750, 500, 250, 150, 100, 75, 50, 40, 30, 20, 10]#distance to next regular puls
    extras=[2000, 1000, 500, 250, 150, 100, 50, 40, 30, 20, 15]#distance to next regular puls
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
    
    plt.plot(extrasFinal, recoveryCycle, marker=".", label="Simulation")
    realDataCMi=[0.324487317, 0.249044214, -0.176721166, 0.452664143, 1.463325936, 1.798056359, 3.110076518, 2.838785715, 3.442295076, 3.763976836,  3.681834103]
    plt.plot(extras, realDataCMi, marker=".", label="CMi")
    realDataCM=[0.240241491, 0.252928694, 0.296115416, 0.546602641, 1.450783851, 2.515537532, 3.780695111, 4.167464007, 4.516383778, 4.916885365, 5.386598886]
    plt.plot(extras, realDataCM, marker=".", label="CM")
    plt.legend()
    plt.xlabel("interspike interval")
    plt.ylabel("slowing/speeding")
    

def plotRealData(filename):
    
    data2, CM, CMi, VHT = dataProcessing.getRealData(filename)

    mean=data2[CMi].mean(axis=1)
    std=data2[CMi].sem(axis=1)

    plt.errorbar(range(len(mean)),mean, std, linestyle='None', marker='.', label="CMi")
    mean=data2[CM].mean(axis=1)
    std=data2[CM].sem(axis=1)
    plt.errorbar(range(len(mean)),mean, std, linestyle='None', marker='.', label="CM")
    mean=data2[VHT].mean(axis=1)
    std=data2[VHT].sem(axis=1)
    plt.errorbar(range(len(mean)),mean, std, linestyle='None', marker='.', label="VHT")
    plt.legend()

#for ELID und Serra
def plotRealSim(filenameReal, data_aps, data_stim, fiberType, norm=False):
    
    data2, CM, CMi, VHT = dataProcessing.getRealData(filenameReal)
    l = dataProcessing.calculateLatency(data_aps, data_stim, norm)
    
    plt.figure()
    
    if fiberType=="CM":
        mean=data2[CM].mean(axis=1)
        std=data2[CM].sem(axis=1)
        plt.plot(l[:,0],mean[0:len(l[:,0])],label="CM", color="red", marker='.')
    elif fiberType=="CMi":
        mean=data2[CMi].mean(axis=1)
        std=data2[CMi].sem(axis=1)
        plt.plot(l[:,0],mean[0:len(l[:,0])],label="CMi", color="red", marker='.')
    elif fiberType=="VHT":
        mean=data2[VHT].mean(axis=1)
        std=data2[VHT].sem(axis=1)
        plt.plot(l[:,0],mean[0:len(l[:,0])],label="VHT", color="red", marker='.')
           
    if len(l)>0:
        plt.plot(l[:,0],l[:,1], label="Simulation", marker='.', zorder=10, color="darkblue")
        plt.xlabel('time (s)')
        if norm:
            plt.ylabel('latency (%)')
        else:
            plt.ylabel('latency (ms)')
        plt.title('Latency')
        plt.legend()
        plt.show()
    
def plotRealSimODP(filenameReal, pathSim, scalingFactor, norm, particle, fibertype="CMi", title="Recovery Cycle"):
    fig, ax=plt.subplots(figsize=(15,5))

    data2, CM, CMi, VHT=dataProcessing.getRealDataODP(filenameReal)
    if fibertype=="CMi":
        mean=data2.iloc[CMi].mean()
        sem=data2.iloc[CMi].sem()
    elif fibertype=="CM":
        mean=data2.iloc[CM].mean()
        sem=data2.iloc[CM].sem()
    elif fibertype=="VHT":
        mean=data2.iloc[VHT].mean()
        sem=data2.iloc[VHT].sem()
    plt.errorbar(range(len(mean)),mean.iloc[::-1], sem.iloc[::-1], linestyle='None', marker='.', label="in-vivo", color="red", ecolor="red")
    
    mean=data2.iloc[CM].mean()
    sem=data2.iloc[CM].sem()
    #plt.errorbar(range(len(mean)),mean.iloc[::-1], sem.iloc[::-1], linestyle='None', marker='.', label="CM", color="blue", ecolor="lightblue")

    data_aps=dataProcessing.getData(path=pathSim, filetype="spikes", prot='DP', scalingFactor=scalingFactor, nr=particle)
    data_stim=pd.read_csv(str(pathSim)+"/stim_ProtDP.csv")

    l=dataProcessing.calculateLatency(data_aps, data_stim, norm=norm)

    recoveryCycle=[]
    for i in range(19,len(l)-2,6):
        recoveryCycle.append(l[i+2][1]-l[i][1])
        #print("Shift")
        #print(i)
        #print(i+2)


    plt.scatter(range(len(recoveryCycle)), recoveryCycle[::-1], label="in-silico", color="darkblue")
    ax.set_xticks(range(11))
    labels=["2000", "1000", "500", "250", "150", "100", "50", "40", "30", "20", "15"]
    ax.set_xticklabels(labels[::-1])
    plt.xlabel("interstimulus interval in ms")
    if norm:
        plt.ylabel('slowing (%)')
    else:
        plt.ylabel('slowing (ms)')
    plt.title(title)
    plt.legend()
    plt.savefig("Figures/Optimization-RecoveryCycle", bbox_inches='tight', dpi=300)
                
def plotRealSimFF(pathSim, scalingFactor, freq, particleNr="best"):
    fig, ax=plt.subplots(figsize=(15,5))

    #Werte von Vivien
    if freq==10:
        #10HZ,CMi, 2-4-8
        mean=[0.92, 2.26, 4.52]
        sem=[0.56, 0.97, 1.56]
    elif freq==20:
        mean=[1.42,3.07,5.48]
        sem=[0.2, 0.16, 1.35]
    elif freq==50:
        mean=[1.27,2.47,4.57]
        sem=[0.27,0.64,1.78]
    elif freq==100:
        mean=[0.84, 0.91,0.99]
        sem=[0.64,0.69,0.57]
        
    plt.errorbar(range(len(mean)),mean, sem, linestyle='None', marker='.', label="CMi", color="black", ecolor="grey")

    data_aps=dataProcessing.getData(path=pathSim, filetype="spikes", prot='FF_'+str(freq)+'Hz_ADS', scalingFactor=scalingFactor, particleNr=particleNr)
    data_stim=dataProcessing.getData(path=pathSim, filetype="stim", prot='FF_'+str(freq)+'Hz_ADS', scalingFactor=scalingFactor)
    #print(data_aps.iloc[54])
    #print(data_stim.iloc[54])
    l=dataProcessing.calculateLatency(data_aps, data_stim, norm=True)
    #print(l)
    
    a=np.where(l[:,0] == 80)
    b=np.where(l[:,0] == 84)
    
    c=np.where(l[:,0] == 120)
    d=np.where(l[:,0] == 124)
    
    e=np.where(l[:,0] == 160)
    f=np.where(l[:,0] == 164)

    if np.isnan(l[f][0,1]):
        print('no AP')
        f=np.where(l[:,0] == 168)

    ff=[l[b][0,1]-l[a][0,1], l[d][0,1]-l[c][0,1], l[f][0,1]-l[e][0,1]]
    print(ff)

    plt.scatter(range(len(ff)), ff, label="Simulation", color="red")
    ax.set_xticks(range(len(ff)))
    labels=[2,4,8]
    ax.set_xticklabels(labels)
    plt.xlabel("number of extra pulses")
    plt.ylabel("slowing in %")
    plt.title("Following Frequencies - "+str(freq)+"Hz, length "+str(scalingFactor*10)+"cm")
    plt.legend()
        
def plotGating(data, title=""):
    data["Nav1.7_m"].plot()
    data["Nav1.7_h"].plot()
    data["Nav1.7_s"].plot()
    plt.legend()
    plt.title(title)
    plt.show()

    data["Nav1.8_m"].plot()
    data["Nav1.8_h"].plot()
    data["Nav1.8_s"].plot()
    data["Nav1.8_u"].plot()
    plt.legend()
    plt.title(title)
    plt.show()

    data["Nav1.9_m"].plot()
    data["Nav1.9_h"].plot()
    data["Nav1.9_s"].plot()
    plt.legend()
    plt.title(title)
    plt.show()

    data["Ks_nf"].plot()
    data["Ks_ns"].plot()
    plt.legend()
    plt.title(title)
    plt.show()

    data["Kf_m"].plot()
    data["Kf_h"].plot()
    plt.legend()
    plt.title(title)
    plt.show()

    data["Kdr_n"].plot()
    plt.legend()
    plt.title(title)
    plt.show()

    data["h_nf"].plot()
    data["h_ns"].plot()
    plt.legend()
    plt.title(title)
    plt.show()        
        
        
        
        