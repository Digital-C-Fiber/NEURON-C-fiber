from neuron import h, gui

def recordV(section):
    v = h.Vector()
    v.record(section._ref_v)
    return v

def recordT():
    t = h.Vector()
    t.record(h._ref_t)
    return t

#didn't work
def recordAP(section):
    spTimes = h.Vector()
    apc = h.APCount(section)
    apc.thresh = -10
    apc.record(spTimes)
    return spTimes