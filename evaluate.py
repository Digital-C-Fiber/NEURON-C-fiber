import main
import dataProcessing

#run model and evaluate goodness of fit

def runParticle(gPump=0.0047891, gNav17=0.10664, gNav18=0.24271, gNav19=9.4779e-05, 
        gKs=0.0069733, gKf=0.012756, gH=0.0025377, gKdr=0.018002, gKna=0.00042):
    prot=-1
    scaling=0.1
    #run model
    #attention: pump has to be negative!
    #run one pulse to check if this gives one AP
    main.run(prot=1, scalingFactor=scaling, gPump=-gPump, gNav17=gNav17, gNav18=gNav18, gNav19=gNav19, gKs=gKs, gKf=gKf, gH=gH, gKdr=gKdr, gKna=gKna)
    dataPulse = dataProcessing.getData(prot=1, filetype="spikes", gPump=-gPump, gNav17=gNav17, gNav18=gNav18, gNav19=gNav19, gKs=gKs, gKf=gKf, gH=gH, gKdr=gKdr, gKna=gKna)
    #if there is only one AP, proceed 
    if dataPulse is not None and dataPulse.shape[0]==1:
        #run test protocol
        main.run(prot=prot, scalingFactor=scaling, gPump=-gPump, gNav17=gNav17, gNav18=gNav18, gNav19=gNav19, gKs=gKs, gKf=gKf, gH=gH, gKdr=gKdr, gKna=gKna)

        return 1
    else:
        return 0

#evaluate goodness of fit of data for specific protocol
def evalParticle(args):
    gPump, gNav17, gNav18, gNav19, gKs, gKf, gH, gKdr, gKna=args
    
    prot=-1 
    #get data
    data_aps = dataProcessing.getData(prot=prot, filetype="spikes", gPump=-gPump, gNav17=gNav17, gNav18=gNav18, gNav19=gNav19, gKs=gKs, gKf=gKf, gH=gH, gKdr=gKdr, gKna=gKna)
    data_stim = dataProcessing.getData(prot=prot, filetype="stim", gPump=-gPump, gNav17=gNav17, gNav18=gNav18, gNav19=gNav19, gKs=gKs, gKf=gKf, gH=gH, gKdr=gKdr, gKna=gKna)
    
    #check if dataframe exists
    if data_aps is not None:
    
        if prot == -1:#test for particle swarm
            #beginning of ELID (20 pulses 1/8Hz), compare ADS with data
            #ADS real: CMi=1.4763, CM=0.0500, difference=1.4263
            #ADS has to go down by 1.4263 to simulate CM
            #ADS sim: CMi=1.4651030783916497
            latency=dataProcessing.calculateLatency(data_aps, data_stim)
            print("Latency: "+str(latency))
            '''
            adsSimCMi=1.4651030783916497
            #y=adsSimCMi-latency[19]
            y=adsSimCMi-latency[4]
            print("to minimize: "+str(y-1.4263))
            return abs(y-1.4263)
            '''
            return abs(latency[4])
        elif prot == 1:
            # compare number of APs: there should only be one
            print("AP number:"+str(data_aps.shape[0]))
            return data_aps.shape[0]
        '''    
        elif prot == 11:#ELID with recovery
            #compare ADS etc

        elif prot == 15:#2HZ with recovery
            #TODO
        '''
    return float('inf')

    