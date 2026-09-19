# plot comprison for ADS and recovery cycle

from neuron import h
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from ast import literal_eval


def getLatency(data_aps, data_stim):
    l = np.zeros((len(data_stim), 2))
    j = 0
    i = 0
    k = 0
    while i < len(data_stim):
        if j < len(data_aps):
            point = data_aps["Axon 3 1"][j]-data_stim["StimTime"][i]
            if point > 0 and point < 400:
                l[k][0] = data_stim["StimTime"][i]/1000
                l[k][1] = point

                j = j+1
                i = i+1
                k = k+1
            elif point < 0:
                j = j+1
            elif point > 400:
                l[k][0] = data_stim["StimTime"][i]/1000
                l[k][1] = -1
                i = i+1
                k = k+1

        else:
            break
    return l


####################### Recovery cycle Comparison (prot16) #######################

# calculate interval-latency change
def getRecoveryCycle(data_aps, data_stim):
    l = getLatency(data_aps, data_stim)
    recoveryCycle = []         # latency change
    extrasFinal = []           # recovery interval
    extras = [
        2000,                       # prot22 comment
        1750,
        1500,
        1250, 1000,
        750, 500, 250,
        150, 100, 75, 50,
        40, 30, 20, 10
    ]
    numRegPulses = 6            # block size: 1 extra + 5 regular pulses
    for i in range(len(extras)):
        l1 = l[20 + i*numRegPulses][1]    # latency of first extra pulse
        l2 = l[21 + i*numRegPulses][1]    # latency of next regular pulse

        if l1 != -1 and l2 != -1:
            recoveryCycle.append(l2 - l1)
            extrasFinal.append(extras[i])

    return np.array(extrasFinal), np.array(recoveryCycle)


def plotRecoveryCompare(original_aps, original_stim,
                        combined_aps, combined_stim):

    x_orig, y_orig = getRecoveryCycle(original_aps, original_stim)
    x_com, y_com = getRecoveryCycle(combined_aps, combined_stim)

    # sort
    order_orig = np.argsort(x_orig)
    order_mod = np.argsort(x_com)

    plt.figure(figsize=(8, 4))
    plt.plot(
        x_orig[order_orig],
        y_orig[order_orig],
        marker="o",
        label="Original"
    )
    plt.plot(
        x_com[order_mod],
        y_com[order_mod],
        marker="o",
        label="3-channel modified",
        color="tab:pink"
    )
    plt.axhline(0, linestyle="--")
    plt.xlim(0, 250)
    plt.ylim(-2, 7)
    plt.xlabel("Interspike interval (ms)")
    plt.ylabel("Latency change (ms)")
    plt.title("Recovery cycle comparison (0.6Hz)")
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":  # data address, need to copy from results
    original_spikes_file = r"C:\Internproj\original\NEURON-C-fiber-main\Results\spikes_Prot24_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"
    original_stim_file = r"C:\Internproj\original\NEURON-C-fiber-main\Results\stim_Prot24_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"
    combined_spikes_file = r"Results\spikes_Prot24_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"
    combined_stim_file = r"Results\stim_Prot24_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"

    original_aps = pd.read_csv(original_spikes_file)
    original_stim = pd.read_csv(original_stim_file)
    combined_aps = pd.read_csv(combined_spikes_file)
    combined_stim = pd.read_csv(combined_stim_file)

    plotRecoveryCompare(
        original_aps,
        original_stim,
        combined_aps,
        combined_stim
    )


####################### Lantency Comparison (prot3,11) #######################
def getLatencyForPlot(data_aps, data_stim):  #
    l = getLatency(data_aps, data_stim)
    #
    valid = l[:, 1] > 0
    latency = l[valid, 1]

    input_number = np.where(valid)[0] + 1

    return input_number, latency

# relative latency


def getRelativeLtency(latency):
    first_latency = latency[0]
    relative_latency = ((latency - first_latency) / first_latency) * 100

    return relative_latency

# comparison


def plotLatencyCompare(original_aps, original_stim,
                       combined_aps, combined_stim):
    # original
    x_orig, latency_orig = getLatencyForPlot(original_aps, original_stim)
    relative_orig = getRelativeLtency(latency_orig)
    # modified
    x_com, latency_com = getLatencyForPlot(combined_aps, combined_stim)
    relative_com = getRelativeLtency(latency_com)

    # A: draw absolute latency
    plt.figure(figsize=(7, 5))
    plt.plot(
        x_orig,
        latency_orig,
        label="Original"
    )
    plt.plot(
        x_com,
        latency_com,
        label="3-channel modified",
        color="tab:pink"
    )
    plt.xlabel("Number of input")
    plt.ylabel("Latency (ms)")
    plt.title("Latency comparison")
    plt.legend()
    plt.grid()

    # B: draw relative latency
    plt.figure(figsize=(7, 5))
    plt.plot(
        x_orig,
        relative_orig,
        label="Original"
    )
    plt.plot(
        x_com,
        relative_com,
        label="3-channel modified",
        color="tab:pink"
    )
    plt.xlabel("Number of input")
    plt.ylabel("Relative latency (%)")
    plt.title("Relative latency comparison")
    plt.legend()
    plt.grid()

    plt.show()


#
# if __name__ == "__main__":  # data address, need to copy from results
#     original_spikes_file = r"C:\Internproj\original\NEURON-C-fiber-main\Results\spikes_Prot11_gPump-0.0047891_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"
#     original_stim_file = r"C:\Internproj\original\NEURON-C-fiber-main\Results\stim_Prot11_gPump-0.0047891_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"
#     combined_spikes_file = r"Results\spikes_Prot11_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"
#     combined_stim_file = r"Results\stim_Prot11_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"

#     original_aps = pd.read_csv(original_spikes_file)
#     original_stim = pd.read_csv(original_stim_file)
#     combined_aps = pd.read_csv(combined_spikes_file)
#     combined_stim = pd.read_csv(combined_stim_file)

#     plotLatencyCompare(
#         original_aps,
#         original_stim,
#         combined_aps,
#         combined_stim
#     )
