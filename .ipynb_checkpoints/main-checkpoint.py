from neuron import h
#from matplotlib import pyplot as plt
import numpy as np
import csv
import time
import math
import os

from defineCell import *
from stimulationProtocols import *
#from saveData import *
#from plot import *

#prot: stimulation protocol number or filename
#gPump: conductance of pump
#gNav17: conductance of Nav 1.7
#gNav18: conductance of Nav 1.8
#dt: step size in time, if set to zero, CVode is activated
#previousStim: sets a pre stimulation before the regular stimulation protocol, if the protocol is loaded from file
def run(prot=1, path="Results/", scalingFactor=1,  dt=0, previousStim=False, tempBranch=32, tempParent=37, 
        gPump=-0.0047891*0.7, gNav17=0.10664, gNav18=0.24271, gNav19=9.4779e-05*20, 
        gKs=0.0069733, gKf=0.012756, gH=0.0025377, gKdr=0.018002*1.2, gKna=0.00042*1.8,vRest=-55,
        sine=False, ampSine=0.1):
    
    #start timer
    tic = time.perf_counter()

    #what do you want to save?
    savePotential=True
    saveCurrents=True
    saveGating=False
    saveSpikes=True
    saveConcentrations=False
    saveStimulation=True

    #create folder
    if not os.path.exists(path):
        os.makedirs(path)
    
    #define morphology as in Tigerholm
    axon=[0,0,0,0,0,0]
    
    axon[0] = h.Section(name='extra1')
    axon[0].nseg = round(10*scalingFactor)
    axon[0].diam = 0.25
    axon[0].L = 100*scalingFactor
    
    axon[1] = h.Section(name='branch')
    axon[1].diam = 0.25
    axon[1].L = 20000*scalingFactor
    axon[1].nseg = round(400*scalingFactor)
        
    axon[2] = h.Section(name='branchPoint')
    h.pt3dadd(h.Vector([0,5000*scalingFactor]), h.Vector([0,0]), h.Vector([0,0]), h.Vector([0.25,1]), sec=axon[2])
    axon[2].nseg = round(100*scalingFactor)
        
    axon[3] = h.Section(name='parent')
    axon[3].diam = 1
    axon[3].L = 100000*scalingFactor
    axon[3].nseg = round(10*200*scalingFactor)
        
    axon[4] = h.Section(name='extra2')
    axon[4].nseg = round(10*scalingFactor)
    axon[4].diam = 1
    axon[4].L = 100*scalingFactor
    
    axon[5] = h.Section(name='extra3')
    axon[5].nseg = round(10*scalingFactor)
    axon[5].diam = 1
    axon[5].L = 100*scalingFactor
    
    for i in range(6):
        #axon[i].Ra = 35.5
        axon[i].Ra = 100
        axon[i].cm = 1
    
    #connect parts
    axon[1].connect(axon[0])
    axon[2].connect(axon[1])
    axon[3].connect(axon[2])
    axon[4].connect(axon[3])
    axon[5].connect(axon[4])
    
    '''
    #print all axon information to check if everything is right
    for i in range(6):
        print(axon[i].psection())
    '''
    
    '''
    h.topology()
    ps = h.PlotShape(False)  # False tells h.PlotShape not to use NEURON's gui
    ps.plot(plt)
    plt.show(0)
    '''
    
    for i in range(6):
        condFactor=1
        if i==0 or i==5:
            condFactor=1e-5
        insertChannels(axon[i], condFactor, gPump, gNav17, gNav18, gNav19, gKs, gKf, gH, gKdr, gKna)
    
    if dt==0:
        #variable time step integration method
        cvode = h.CVode()
        cvode.active(1)
        print("CVode active")
    else:
        h.dt=dt

    if sine:
        stim, delay, vec = setStimulationSine(axon[0], prot, ampSine)
        print("Sine")
    else:
        stim, delay, vec = setStimulationProtocol(axon[0], prot, previousStim)

    if saveSpikes:
        spTimes = h.Vector()
        apc = h.APCount(axon[1](0))
        apc.thresh = -10
        apc.record(spTimes)
    
        spTimes2 = h.Vector() 
        apc2 = h.APCount(axon[1](0.25))
        apc2.thresh = -10
        apc2.record(spTimes2)
        
        spTimes3 = h.Vector() 
        apc3 = h.APCount(axon[1](0.5))
        apc3.thresh = -10
        apc3.record(spTimes3)
        
        spTimes4 = h.Vector() 
        apc4 = h.APCount(axon[1](0.75))
        apc4.thresh = -10
        apc4.record(spTimes4)
        
        spTimes5 = h.Vector() 
        apc5 = h.APCount(axon[1](1))
        apc5.thresh = -10
        apc5.record(spTimes5)
        
        spTimes6 = h.Vector() 
        apc6 = h.APCount(axon[3](0))
        apc6.thresh = -10
        apc6.record(spTimes6)
        
        spTimes7 = h.Vector() 
        apc7 = h.APCount(axon[3](0.25))
        apc7.thresh = -10
        apc7.record(spTimes7)
        
        spTimes8 = h.Vector() 
        apc8 = h.APCount(axon[3](0.5))
        apc8.thresh = -10
        apc8.record(spTimes8)
        
        spTimes9 = h.Vector() 
        apc9 = h.APCount(axon[3](0.75))
        apc9.thresh = -10
        apc9.record(spTimes9)
        
        spTimes10 = h.Vector() 
        apc10 = h.APCount(axon[3](1))
        apc10.thresh = -10
        apc10.record(spTimes10)

    #record time
    t = h.Vector().record(h._ref_t)
    
    if savePotential:    
        potential = h.Vector().record(axon[3](1)._ref_v)
        
    if saveCurrents:
        #record currents
        nav17 = h.Vector()
        nav17.record(axon[3](1)._ref_ina_nattxs)
        nav18 = h.Vector()
        nav18.record(axon[3](1)._ref_ina_nav1p8)
        nav19 = h.Vector()
        nav19.record(axon[3](1)._ref_ina_nav1p9)
        ks = h.Vector()
        ks.record(axon[3](1)._ref_ik_ks)
        kf = h.Vector()
        kf.record(axon[3](1)._ref_ik_kf)
        kdr = h.Vector()
        kdr.record(axon[3](1)._ref_ik_kdr)
        h_na = h.Vector()
        h_na.record(axon[3](1)._ref_ina_h)
        h_k = h.Vector()
        h_k.record(axon[3](1)._ref_ik_h)
        pump_na = h.Vector()
        pump_na.record(axon[3](1)._ref_ina_nakpump)
        pump_k = h.Vector()
        pump_k.record(axon[3](1)._ref_ik_nakpump)
        leak_na = h.Vector()
        leak_na.record(axon[3](1)._ref_ina_leak)
        leak_k = h.Vector()
        leak_k.record(axon[3](1)._ref_ik_leak)
        kna = h.Vector()
        kna.record(axon[3](1)._ref_ik_kna)
    
    if saveGating:
        #record gating variables
        nav17_m = h.Vector()
        nav17_m.record(axon[3](1)._ref_m_nattxs)
        nav17_h = h.Vector()
        nav17_h.record(axon[3](1)._ref_h_nattxs)
        nav17_s = h.Vector()
        nav17_s.record(axon[3](1)._ref_s_nattxs)

        nav18_m = h.Vector()
        nav18_m.record(axon[3](1)._ref_m_nav1p8)
        nav18_h = h.Vector()
        nav18_h.record(axon[3](1)._ref_h_nav1p8)
        nav18_s = h.Vector()
        nav18_s.record(axon[3](1)._ref_s_nav1p8)
        nav18_u = h.Vector()
        nav18_u.record(axon[3](1)._ref_u_nav1p8)

        nav19_m = h.Vector()
        nav19_m.record(axon[3](1)._ref_m_nav1p9)
        nav19_h = h.Vector()
        nav19_h.record(axon[3](1)._ref_h_nav1p9)
        nav19_s = h.Vector()
        nav19_s.record(axon[3](1)._ref_s_nav1p9)

        ks_nf = h.Vector()
        ks_nf.record(axon[3](1)._ref_nf_ks)
        ks_ns = h.Vector()
        ks_ns.record(axon[3](1)._ref_ns_ks)

        kf_m = h.Vector()
        kf_m.record(axon[3](1)._ref_m_kf)
        kf_h = h.Vector()
        kf_h.record(axon[3](1)._ref_h_kf)

        kdr_n = h.Vector()
        kdr_n.record(axon[3](1)._ref_n_kdr)

        h_nf = h.Vector()
        h_nf.record(axon[3](1)._ref_nf_h)
        h_ns = h.Vector()
        h_ns.record(axon[3](1)._ref_ns_h)
        ''' 
        #Branch
        t1 = h.Vector().record(h._ref_t)
        potential1 = h.Vector().record(axon[1](1)._ref_v)
        
        #record currents
        nav17_1 = h.Vector()
        nav17_1.record(axon[1](1)._ref_ina_nattxs)
        nav18_1 = h.Vector()
        nav18_1.record(axon[1](1)._ref_ina_nav1p8)
        nav19_1 = h.Vector()
        nav19_1.record(axon[1](1)._ref_ina_nav1p9)
        ks_1 = h.Vector()
        ks_1.record(axon[1](1)._ref_ik_ks)
        kf_1 = h.Vector()
        kf_1.record(axon[1](1)._ref_ik_kf)
        kdr_1 = h.Vector()
        kdr_1.record(axon[1](1)._ref_ik_kdr)
        h_na_1 = h.Vector()
        h_na_1.record(axon[1](1)._ref_ina_h)
        h_k_1 = h.Vector()
        h_k_1.record(axon[1](1)._ref_ik_h)
        pump_na_1 = h.Vector()
        pump_na_1.record(axon[1](1)._ref_ina_nakpump)
        pump_k_1 = h.Vector()
        pump_k_1.record(axon[1](1)._ref_ik_nakpump)
        leak_na_1 = h.Vector()
        leak_na_1.record(axon[1](1)._ref_ina_leak)
        leak_k_1 = h.Vector()
        leak_k_1.record(axon[1](1)._ref_ik_leak)
        kna_1 = h.Vector()
        kna_1.record(axon[1](1)._ref_ik_kna)

        #record gating variables
        nav17_m_1 = h.Vector()
        nav17_m_1.record(axon[1](1)._ref_m_nattxs)
        nav17_h_1 = h.Vector()
        nav17_h_1.record(axon[1](1)._ref_h_nattxs)
        nav17_s_1 = h.Vector()
        nav17_s_1.record(axon[1](1)._ref_s_nattxs)

        nav18_m_1 = h.Vector()
        nav18_m_1.record(axon[1](1)._ref_m_nav1p8)
        nav18_h_1 = h.Vector()
        nav18_h_1.record(axon[1](1)._ref_h_nav1p8)
        nav18_s_1 = h.Vector()
        nav18_s_1.record(axon[1](1)._ref_s_nav1p8)
        nav18_u_1 = h.Vector()
        nav18_u_1.record(axon[1](1)._ref_u_nav1p8)

        nav19_m_1 = h.Vector()
        nav19_m_1.record(axon[1](1)._ref_m_nav1p9)
        nav19_h_1 = h.Vector()
        nav19_h_1.record(axon[1](1)._ref_h_nav1p9)
        nav19_s_1 = h.Vector()
        nav19_s_1.record(axon[1](1)._ref_s_nav1p9)

        ks_nf_1 = h.Vector()
        ks_nf_1.record(axon[1](1)._ref_nf_ks)
        ks_ns_1 = h.Vector()
        ks_ns_1.record(axon[1](1)._ref_ns_ks)

        kf_m_1 = h.Vector()
        kf_m_1.record(axon[1](1)._ref_m_kf)
        kf_h_1 = h.Vector()
        kf_h_1.record(axon[1](1)._ref_h_kf)

        kdr_n_1 = h.Vector()
        kdr_n_1.record(axon[1](1)._ref_n_kdr)

        h_nf_1 = h.Vector()
        h_nf_1.record(axon[1](1)._ref_nf_h)
        h_ns_1 = h.Vector()
        h_ns_1.record(axon[1](1)._ref_ns_h)
        '''
    if saveConcentrations:
        save_nai = h.Vector().record(axon[3](1)._ref_nai)
        save_ki = h.Vector().record(axon[3](1)._ref_ki)
        save_nao = h.Vector().record(axon[3](1)._ref_nao)
        save_ko = h.Vector().record(axon[3](1)._ref_ko)
        save_ena = h.Vector().record(axon[3](1)._ref_ena)
        save_ek = h.Vector().record(axon[3](1)._ref_ek)

    #simulation
    Vrest=vRest
    h.finitialize(Vrest)

    tempCelsius1 = tempBranch
    tempCelsius2 = tempParent
    setTemp(axon[0], tempCelsius1)    
    setTemp(axon[1], tempCelsius1)
    setTemp(axon[2], (tempCelsius1+tempCelsius2)/2)
    setTemp(axon[3], tempCelsius2)
    setTemp(axon[4], tempCelsius2)
    setTemp(axon[5], tempCelsius2)
    

    h.fcurrent()
    
    for i in range(6):
        balance(axon[i], Vrest)
    
    #create filename
    if isinstance(prot, str):
        prot = prot.split("/", 2)
        prot = prot[2]
    
    #filename can't be too long, full path can't be more than 255 characters
    '''
    fileSuffix=('_Prot'+str(prot)+'_scalingFactor'+str(scalingFactor)
                +'_tempBranch'+str(tempBranch)+'_tempParent'+str(tempParent)
                +'_gPump'+str(gPump)+'_gNav17'+str(gNav17)+'_gNav18'+str(gNav18)+'_gNav19'+str(gNav19)
                +'_gKs'+str(gKs)+'_gKf'+str(gKf)+'_gH'+str(gH)+'_gKdr'+str(gKdr)+'_gKna'+str(gKna)+'_vRest'+str(vRest)+'.csv')
    '''
    '''
    fileSuffix=('_Prot'+str(prot)
                +'_gPump'+str(round(gPump,5))+'_gNav17'+str(gNav17)+'_gNav18'+str(gNav18)+'_gNav19'+str(gNav19)
                +'_gKs'+str(gKs)+'_gKf'+str(gKf)+'_gH'+str(gH)+'_gKdr'+str(gKdr)+'_gKna'+str(gKna)+'_vRest'+str(vRest)
                +'.csv')
    '''
    fileSuffix=('_Prot'+str(prot)
                +'_gPump'+str(round(gPump,4))+'_gNav17'+str(round(gNav17,4))+'_gNav18'+str(round(gNav18,4))+'_gNav19'+str(round(gNav19,4))
                +'_gKs'+str(round(gKs,4))+'_gKf'+str(round(gKf,4))+'_gH'+str(round(gH,4))+'_gKdr'+str(round(gKdr,4))+'_gKna'+str(round(gKna,4))+'_vRest'+str(round(vRest,4))
                +'.csv')

    if savePotential:
        filename = path+'potential'+fileSuffix
        
        #creates file, deletes content, if file already exists
        with open(filename,'w', newline='') as f:
            csv.writer(f).writerow(["Time", "Axon 1 0", "Axon 1 0.25", "Axon 1 0.5", "Axon 1 0.75", "Axon 1 1", "Axon 3 0", "Axon 3 0.25", "Axon 3 0.5","Axon 3 0.75","Axon 3 1"])

    if saveSpikes:
        fileSpikes = path+'spikes'+fileSuffix
        with open(fileSpikes,'w', newline='') as f:
            csv.writer(f).writerow(["Axon 1 0", "Axon 1 0.25", "Axon 1 0.5", "Axon 1 0.75", "Axon 1 1", "Axon 3 0", "Axon 3 0.25", "Axon 3 0.5","Axon 3 0.75","Axon 3 1"])

    if saveStimulation:
        #Stimulation times
        fileStim = path+'stim'+fileSuffix
        with open(fileStim,'w', newline='') as f:
            csv.writer(f).writerow(["StimTime"])
            
        for stimTime in vec:
            with open(fileStim,'a', newline='') as f:
                csv.writer(f).writerow([stimTime])
    
    #start simulation
    tstop = delay
    #h.continuerun(tstop)
    i=0
    while(h.t<tstop):
        #save data
        if savePotential:
            with open(filename,'a', newline='') as f:
                csv.writer(f).writerow([h.t, axon[1](0).v, axon[1](0.25).v, axon[1](0.5).v, axon[1](0.75).v, axon[1](1).v, axon[3](0).v, axon[3](0.25).v, axon[3](0.5).v, axon[3](0.75).v, axon[3](1).v])

        if saveSpikes:
            if i<len(spTimes) and i<len(spTimes2) and i<len(spTimes3) and i<len(spTimes4) and i<len(spTimes5) and i<len(spTimes6) and i<len(spTimes7) and i<len(spTimes8) and i<len(spTimes9) and i<len(spTimes10):
                print("Percentage done: "+ str(h.t/tstop*100))
                #print("Time: "+str(h.t))
                #print("APs:")
                #print("Axon 1 0.5: " + str(spTimes3[i]))
                #print("Axon 3 0.5: " + str(spTimes8[i]))
                
                with open(fileSpikes,'a', newline='') as f:
                    
                    csv.writer(f).writerow([spTimes[i], spTimes2[i],spTimes3[i], spTimes4[i], spTimes5[i], spTimes6[i], 
                                            spTimes7[i], spTimes8[i],spTimes9[i], spTimes10[i]])
                    
                i=i+1
            
        #step 
        h.fadvance()
    
    if saveCurrents:
        filename = str(path)+'/cur'+fileSuffix    
        with open(filename,'w', newline='') as f:
            csv.writer(f).writerow(["Time","Nav1.7", "Nav1.8","Nav1.9", "Ks", "Kf", "Kdr", "h_na", "h_k", "pump_na", "pump_k", "leak_na", "leak_k", "kna"])
            csv.writer(f).writerows(zip(*(t, nav17, nav18, nav19, ks, kf, kdr, h_na, h_k, pump_na, pump_k, leak_na, leak_k, kna)))

    if saveGating:
        filename = str(path)+'/gat'+fileSuffix    
        with open(filename,'w', newline='') as f:
            csv.writer(f).writerow(["Time","Nav1.7_m","Nav1.7_h","Nav1.7_s", "Nav1.8_m","Nav1.8_h","Nav1.8_s","Nav1.8_u", "Nav1.9_m", "Nav1.9_h", "Nav1.9_s", "Ks_nf", "Ks_ns", "Kf_m", "Kf_h", "Kdr_n", "h_nf", "h_ns"])
            csv.writer(f).writerows(zip(*(t, nav17_m, nav17_h, nav17_s, nav18_m, nav18_h, nav18_s, nav18_u, nav19_m, nav19_h, nav19_s, ks_nf, ks_ns, kf_m, kf_h, kdr_n, h_nf, h_ns)))
        
    if saveConcentrations:
        filenameConc = str(path)+'/concentrations'+fileSuffix
        #creates file, deletes content, if file already exists
        with open(filenameConc,'w', newline='') as f:
            csv.writer(f).writerow(["Time", "Nai", "Ki", "Nao", "Ko", "ena", "ek"])
            csv.writer(f).writerows(zip(*(t, save_nai, save_ki, save_nao, save_ko, save_ena, save_ek)))

    
    toc = time.perf_counter()
    print(f"Simulation time: {(toc - tic)/60:0.4f} min")


