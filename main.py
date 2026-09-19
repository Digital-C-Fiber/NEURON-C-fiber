from neuron import h
# from matplotlib import pyplot as plt
import numpy as np
import csv
import time
import math
import os
from os import path

from defineCell import *
from stimulationProtocols import *
# from saveData import *
# from plot import *

# prot: stimulation protocol number or filename
# gPump: conductance of pump
# gNav13: conductance of Nav 1.3
# gNav17: conductance of Nav 1.7
# gNav18: conductance of Nav 1.8
# dt: step size in time, if set to zero, CVode is activated
# previousStim: sets a pre stimulation before the regular stimulation protocol, if the protocol is loaded from file


def run(prot=1, scalingFactor=1,  dt=0, previousStim=False, tempBranch=32, tempParent=37,
        gPump=-0.0047891, gNav13=0.010664, gNav17=0.10664, gNav18=0.24271, gNav19=9.4779e-05,
        gKs=0.0069733, gKf=0.012756, gH=0.0025377, gKdr=0.018002, gKna=0.00042, vRest=-55,
        sine=False, ampSine=0.1):

    # start timer
    tic = time.perf_counter()

    # define morphology as in Tigerholm
    axon = [0, 0, 0, 0, 0, 0]

    axon[0] = h.Section(name='extra1')
    axon[0].nseg = round(10*scalingFactor)
    axon[0].diam = 0.25
    axon[0].L = 100*scalingFactor

    axon[1] = h.Section(name='branch')
    axon[1].diam = 0.25
    axon[1].L = 20000*scalingFactor
    axon[1].nseg = round(400*scalingFactor)

    axon[2] = h.Section(name='branchPoint')
    h.pt3dadd(h.Vector([0, 5000*scalingFactor]), h.Vector([0, 0]),
              h.Vector([0, 0]), h.Vector([0.25, 1]), sec=axon[2])
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
        axon[i].Ra = 35.5
        axon[i].cm = 1

    # connect parts
    axon[1].connect(axon[0])
    axon[2].connect(axon[1])
    axon[3].connect(axon[2])
    axon[4].connect(axon[3])
    axon[5].connect(axon[4])

    # '''
    # print all axon information to check if everything is right
    # for i in range(6):
    #    print(axon[i].psection())
    # '''

    # '''
    # h.topology()
    # ps = h.PlotShape(False)  # False tells h.PlotShape not to use NEURON's gui
    # ps.plot(plt)
    # plt.show(0)
    # '''

    for i in range(6):
        condFactor = 1
        if i == 0 or i == 5:
            condFactor = 1e-5
        insertChannels(axon[i], condFactor, gPump, gNav13, gNav17,
                       gNav18, gNav19, gKs, gKf, gH, gKdr, gKna)

    if dt == 0:
        # variable time step integration method
        cvode = h.CVode()
        cvode.active(1)
        print("CVode active")
    else:
        h.dt = dt

    if sine:
        stim, delay, vec = setStimulationSine(axon[0], prot, ampSine)
        print("Stimulation: Sine Wave")
    else:
        stim, delay, vec = setStimulationProtocol(axon[0], prot, previousStim)
        print("Stimulation: Square Pulse")

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

    # Recording
    ############################################################
    rec_seg = axon[3](0.5)

    t_total = h.Vector().record(h._ref_t)
    v_total = h.Vector().record(rec_seg._ref_v)
    ina_total = h.Vector().record(rec_seg._ref_ina)

    ### 1.7 ###
    t_nav17 = h.Vector().record(h._ref_t)
    v_nav17 = h.Vector().record(rec_seg._ref_v)
    ina17 = h.Vector().record(rec_seg._ref_ina_nattxs)

    m17 = h.Vector().record(rec_seg.nattxs._ref_m)
    s17 = h.Vector().record(rec_seg.nattxs._ref_s)

    h_eff17 = h.Vector().record(rec_seg.nattxs._ref_h_eff)
    h0_17 = h.Vector().record(rec_seg.nattxs_h._ref_h0)
    h1_17 = h.Vector().record(rec_seg.nattxs_h._ref_h1)
    hw1_17 = h.Vector().record(rec_seg.nattxs_h._ref_hw1)
    hw2_17 = h.Vector().record(rec_seg.nattxs_h._ref_hw2)
    hw3_17 = h.Vector().record(rec_seg.nattxs_h._ref_hw3)
    hw4_17 = h.Vector().record(rec_seg.nattxs_h._ref_hw4)
    hw5_17 = h.Vector().record(rec_seg.nattxs_h._ref_hw5)

    ### 1.8 ###
    t_nav18 = h.Vector().record(h._ref_t)
    v_nav18 = h.Vector().record(rec_seg._ref_v)
    ina18 = h.Vector().record(rec_seg._ref_ina_nav1p8)

    m18 = h.Vector().record(rec_seg.nav1p8._ref_m)
    s18 = h.Vector().record(rec_seg.nav1p8._ref_s)
    u18 = h.Vector().record(rec_seg.nav1p8._ref_u)

    h_eff18 = h.Vector().record(rec_seg.nav1p8._ref_h_eff)
    h0_18 = h.Vector().record(rec_seg.nav1p8_h._ref_h0)
    h1_18 = h.Vector().record(rec_seg.nav1p8_h._ref_h1)
    hw1_18 = h.Vector().record(rec_seg.nav1p8_h._ref_hw1)
    hw2_18 = h.Vector().record(rec_seg.nav1p8_h._ref_hw2)
    hw3_18 = h.Vector().record(rec_seg.nav1p8_h._ref_hw3)
    hw4_18 = h.Vector().record(rec_seg.nav1p8_h._ref_hw4)
    hw5_18 = h.Vector().record(rec_seg.nav1p8_h._ref_hw5)

    ### 1.3 ###
    t_nav13 = h.Vector().record(h._ref_t)
    v_nav13 = h.Vector().record(rec_seg._ref_v)
    ina13 = h.Vector().record(rec_seg._ref_ina_nav1p3)

    m13 = h.Vector().record(rec_seg.nav1p3._ref_m)

    h_eff13 = h.Vector().record(rec_seg.nav1p3._ref_h_eff)
    h0_13 = h.Vector().record(rec_seg.nav1p3_h._ref_h0)
    h1_13 = h.Vector().record(rec_seg.nav1p3_h._ref_h1)
    hw1_13 = h.Vector().record(rec_seg.nav1p3_h._ref_hw1)
    hw2_13 = h.Vector().record(rec_seg.nav1p3_h._ref_hw2)
    hw3_13 = h.Vector().record(rec_seg.nav1p3_h._ref_hw3)
    hw4_13 = h.Vector().record(rec_seg.nav1p3_h._ref_hw4)
    hw5_13 = h.Vector().record(rec_seg.nav1p3_h._ref_hw5)

    # simulation
    Vrest = vRest

    tempCelsius1 = tempBranch
    tempCelsius2 = tempParent
    setTemp(axon[0], tempCelsius1)
    setTemp(axon[1], tempCelsius1)
    setTemp(axon[2], (tempCelsius1+tempCelsius2)/2)
    setTemp(axon[3], tempCelsius2)
    setTemp(axon[4], tempCelsius2)
    setTemp(axon[5], tempCelsius2)

    h.finitialize(Vrest)
    h.fcurrent()

    for i in range(6):
        balance(axon[i], Vrest)

    # create folder
    if not os.path.exists('Results'):
        os.mkdir('Results')

    # create filename
    if isinstance(prot, str) and "/" in prot:
        prot = prot.split("/", 2)
        print(prot)
        prot = prot[2]

    # filename can't be too long, full path can't be more than 255 characters
    # therefore values are rounded!
    # '''
    # fileSuffix=('_Prot'+str(prot)+'_scalingFactor'+str(scalingFactor)
    #            +'_tempBranch'+str(tempBranch)+'_tempParent'+str(tempParent)
    #            +'_gPump'+str(gPump)+'_gNav17'+str(gNav17)+'_gNav18'+str(gNav18)+'_gNav19'+str(gNav19)
    #            +'_gKs'+str(gKs)+'_gKf'+str(gKf)+'_gH'+str(gH)+'_gKdr'+str(gKdr)+'_gKna'+str(gKna)+'_vRest'+str(vRest)+'.csv')
    # '''
    fileSuffix = ('_Prot'+str(prot)
                  + '_gPump'+str(round(gPump, 7))
                  + '_gNav13'+str(round(gNav13, 7))
                  + '_gNav17'+str(round(gNav17, 7))
                  + '_gNav18'+str(round(gNav18, 7))
                  + '_gNav19'+str(round(gNav19, 7))
                  + '_gKs'+str(round(gKs, 7))
                  + '_gKf'+str(round(gKf, 7))
                  + '_gH'+str(round(gH, 7))
                  + '_gKdr'+str(round(gKdr, 7))
                  + '_gKna'+str(round(gKna, 7))
                  + '_vRest'+str(vRest)
                  + '_sine'+str(sine)
                  + '_ampSine'+str(ampSine)
                  + '.csv')
    filename = 'Results/potential'+fileSuffix

    fileNav13 = f'Results/combined_nav13.csv'
    fileNav17 = f'Results/combined_nav17.csv'
    fileNav18 = f'Results/combined_nav18.csv'

    fileOverall = 'Results/combined_overall.csv'

    # creates file, deletes content, if file already exists
    with open(filename, 'w', newline='') as f:
        csv.writer(f).writerow(["Time", "Axon 1 0", "Axon 1 0.25", "Axon 1 0.5", "Axon 1 0.75",
                                "Axon 1 1", "Axon 3 0", "Axon 3 0.25", "Axon 3 0.5", "Axon 3 0.75", "Axon 3 1"])

    fileSpikes = 'Results/spikes'+fileSuffix
    with open(fileSpikes, 'w', newline='') as f:
        csv.writer(f).writerow(["Axon 1 0", "Axon 1 0.25", "Axon 1 0.5", "Axon 1 0.75",
                                "Axon 1 1", "Axon 3 0", "Axon 3 0.25", "Axon 3 0.5", "Axon 3 0.75", "Axon 3 1"])

    # Stimulation times
    fileStim = 'Results/stim'+fileSuffix
    with open(fileStim, 'w', newline='') as f:
        csv.writer(f).writerow(["StimTime"])

    for stimTime in vec:
        with open(fileStim, 'a', newline='') as f:
            csv.writer(f).writerow([stimTime])

    # start simulation
    tstop = delay
    # h.continuerun(tstop)
    i = 0
    while (h.t < tstop):
        # save data
        with open(filename, 'a', newline='') as f:
            csv.writer(f).writerow([h.t, axon[1](0).v, axon[1](0.25).v, axon[1](0.5).v, axon[1](0.75).v, axon[1](
                1).v, axon[3](0).v, axon[3](0.25).v, axon[3](0.5).v, axon[3](0.75).v, axon[3](1).v])
        if i < len(spTimes) and i < len(spTimes2) and i < len(spTimes3) and i < len(spTimes4) and i < len(spTimes5) and i < len(spTimes6) and i < len(spTimes7) and i < len(spTimes8) and i < len(spTimes9) and i < len(spTimes10):
            print("Time: "+str(h.t))
            print("AP number:"+str(i+1))
            print("Axon 1 0.5: " + str(spTimes3[i]))
            print("Axon 3 0.5: " + str(spTimes8[i]))
            with open(fileSpikes, 'a', newline='') as f:

                csv.writer(f).writerow([spTimes[i], spTimes2[i], spTimes3[i], spTimes4[i], spTimes5[i], spTimes6[i],
                                        spTimes7[i], spTimes8[i], spTimes9[i], spTimes10[i]])

            i = i+1

        # step
        h.fadvance()

    # Save Nav1.7 recordings after stimulation
    ##########################################
    with open(fileNav17, 'w', newline='') as f:
        writer = csv.writer(f)
        # column names
        writer.writerow([
            "Time", "Vm", "Ina17",
            "m",
            "s",
            "h_eff",
            "h0", "h1",
            "hw1", "hw2", "hw3", "hw4", "hw5"
        ])
        # recorded values
        writer.writerows(zip(
            t_nav17, v_nav17, ina17,
            m17,
            s17,
            h_eff17,
            h0_17, h1_17,
            hw1_17, hw2_17, hw3_17, hw4_17, hw5_17
        ))

    # Save Nav1.8 recordings after stimulation
    ##########################################
    with open(fileNav18, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            "Time", "Vm", "Ina18",
            "m",
            "s",
            "u",
            "h_eff",
            "h0", "h1",
            "hw1", "hw2", "hw3", "hw4", "hw5"
        ])
        writer.writerows(zip(
            t_nav18, v_nav18, ina18,
            m18,
            s18, u18,
            h_eff18,
            h0_18, h1_18,
            hw1_18, hw2_18, hw3_18, hw4_18, hw5_18
        ))

    # Save Nav1.3 recordings after stimulation
    ##########################################
    with open(fileNav13, 'w', newline='') as f:
        writer = csv.writer(f)
        # column names
        writer.writerow([
            "Time", "Vm", "Ina13",
            "m",
            "h_eff",
            "h0", "h1",
            "hw1", "hw2", "hw3", "hw4", "hw5"
        ])
        # recorded values
        writer.writerows(zip(
            t_nav13, v_nav13, ina13,
            m13, h_eff13,
            h0_13, h1_13,
            hw1_13, hw2_13, hw3_13, hw4_13, hw5_13
        ))

    # Save pverall recordings after stimulation
    ###########################################
    with open(fileOverall, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            "Time",
            "Vm",
            "Ina"
        ])
        writer.writerows(zip(
            t_total,
            v_total,
            ina_total
        ))

    # plot
    # l = plotLatency(spTimes10, vec)

    # print(v3)
    # for x in spTimes: print(x)

    toc = time.perf_counter()
    print(f"Simulation time: {(toc - tic)/60:0.4f} min")
    # return t, v_mid


# '''
# def range_assignment(sec, var, start, stop):
#     """linearly assign values between start and stop to each segment.var in section"""
#     import numpy as np
#     for seg, val in zip(sec, np.linspace(start, stop, sec.nseg)):
#         setattr(seg, var, val)
# '''
