from neuron import h

def insertChannels(axon,condFactor, gPump=-0.0047891, gNav17=0.10664, gNav18=0.24271, gNav19=9.4779e-05, gKs=0.0069733, gKf=0.012756, gH=0.0025377, gKdr=0.018002, gKna=0.00042):
    axon.insert('ks')
    axon.insert('kf')
    axon.insert('h')
    axon.insert('nattxs')
    axon.insert('nav1p8')
    axon.insert('nav1p9')
    axon.insert('nakpump')
    axon.insert('kdr')#Tigerholm
    #axon.insert('kdrTiger')#Grill
    axon.insert('kna')
    axon.insert('extrapump')
    axon.insert('leak')
    
    axon.insert('nakdyn')#Tigerholm
    #axon.insert('naoi')#Grill
    #axon.insert('koi')#Grill

    for seg in axon:
        seg.ks.gbar = gKs*condFactor
        seg.kf.gbar = gKf*condFactor
        seg.h.gbar = gH*condFactor
        #seg.nattxs.gbar = 0.10664#Nav1.7
        seg.nattxs.gbar = gNav17*condFactor
        #seg.nav1p8.gbar = 0.24271
        seg.nav1p8.gbar = gNav18*condFactor
        seg.nav1p9.gbar = gNav19*condFactor
        #seg.nakpump.smalla = -0.0047891
        seg.nakpump.smalla = gPump*condFactor
        seg.kdr.gbar = gKdr*condFactor#Tigerholm
        #seg.kdrTiger.gbar = 0.018002*condFactor#Grill
        seg.kna.gbar = gKna*condFactor
        
    h.theta_nakdyn = 0.029#Tigerholm
    #h.theta_naoi = 0.029*condFactor#Grill
    #h.theta_koi = 0.029*condFactor#Grill

def setTemp(axon, tempCelsius):
    for seg in axon:
        seg.ks.celsiusT = tempCelsius
        seg.kf.celsiusT = tempCelsius
        seg.h.celsiusT = tempCelsius
        seg.nattxs.celsiusT = tempCelsius
        seg.nav1p9.celsiusT = tempCelsius
        seg.nakpump.celsiusT = tempCelsius
        seg.kdr.celsiusT = tempCelsius#Tigerholm
        #seg.kdrTiger.celsiusT = tempCelsius#Grill
    h.celsiusT_nav1p8 = tempCelsius

#calculates the maximum conductances of the sodium and potassium leak (balancing) currents to achieve their target rest potential (Grill)
def balance(axon, Vrest):
    inaSum = -(axon.ina_nattxs + axon.ina_nav1p9+axon.ina_nav1p8 + axon.ina_h+axon.ina_nakpump)
    if (inaSum/(Vrest-axon.ena))<0:
        axon.pumpina_extrapump=inaSum
    else:
        axon.gnaleak_leak = inaSum/(Vrest-axon.ena)

    ikSum = -(axon.ik_ks+axon.ik_kf + axon.ik_h + axon.ik_kdr + axon.ik_nakpump + axon.ik_kna)#Tigerholm
    #ikSum = -(axon.ik_ks+axon.ik_kf + axon.ik_h + axon.ik_kdrTiger + axon.ik_nakpump + axon.ik_kna)#Grill
    if (ikSum/(Vrest-axon.ek))<0:
        axon.pumpik_extrapump=ikSum
    else:
        axon.gkleak_leak = ikSum/(Vrest-axon.ek)
       


    
    
    
    
    
    
    