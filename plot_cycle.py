import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def getLatency(data_aps, data_stim):
    l = np.zeros((len(data_stim), 2))
    j = 0
    i = 0
    k = 0

    while i < len(data_stim):
        if j < len(data_aps):
            point = data_aps["Axon 3 1"][j] - data_stim["StimTime"][i]

            if point > 0 and point < 400:
                l[k][0] = data_stim["StimTime"][i] / 1000
                l[k][1] = point
                j += 1
                i += 1
                k += 1

            elif point < 0:
                j += 1

            elif point > 400:
                l[k][0] = data_stim["StimTime"][i] / 1000
                l[k][1] = -1
                i += 1
                k += 1
        else:
            break

    return l


def getRecoveryCycle(data_aps, data_stim, extras, numRegPulses=5):
    l = getLatency(data_aps, data_stim)

    recoveryCycle = []
    extrasFinal = []

    for i in range(len(extras)):
        idx1 = 20 + i * (numRegPulses + 1)   # extra pulse
        idx2 = idx1 + 1                      # first regular pulse after extra

        if idx2 < len(l):
            l1 = l[idx1][1]
            l2 = l[idx2][1]

            if l1 > 0 and l2 > 0:
                recoveryCycle.append(l2 - l1)
                extrasFinal.append(extras[i])

    return np.array(extrasFinal), np.array(recoveryCycle)


def plotRecoveryCompareFreq(aps_04, stim_04, aps_05, stim_05, aps_06, stim_06):
    # same in protocol
    extras_04 = [2000, 1750, 1500, 1250, 1000, 750,
                 500, 250, 150, 100, 75, 50, 40, 30, 20, 10]
    extras_05 = [1750, 1500, 1250, 1000, 750, 500,
                 250, 150, 100, 75, 50, 40, 30, 20, 10]
    extras_06 = [1500, 1250, 1000, 750, 500,
                 250, 150, 100, 75, 50, 40, 30, 20, 10]

    x04, y04 = getRecoveryCycle(aps_04, stim_04, extras_04)
    x05, y05 = getRecoveryCycle(aps_05, stim_05, extras_05)
    x06, y06 = getRecoveryCycle(aps_06, stim_06, extras_06)

    plt.figure(figsize=(8, 4))

    plt.plot(x06, y06, linestyle="-", marker="o",
             markersize=4, color="black",  label="0.6 Hz")
    plt.plot(x05, y05, linestyle="-", marker="o",
             markersize=4, color="red", label="0.5 Hz")
    plt.plot(x04, y04, linestyle="-", marker="o", markersize=4,
             color="darkturquoise",   label="0.4 Hz")

    plt.axhline(0, linestyle="--", color="gray")
    plt.xlabel("Interspike interval (ms)")
    plt.ylabel("Latency change (ms)")
    plt.title("Recovery cycle")
    plt.xlim(0, 250)
    plt.ylim(-1.6, 6.5)
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    # =========================
    # 0.4 Hz
    # =========================
    spikes_04 = r"C:\Internproj\combine\NEURON-C-fiber-main\Results\spikes_Prot25_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"
    stim_04 = r"C:\Internproj\combine\NEURON-C-fiber-main\Results\stim_Prot25_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"

    # =========================
    # 0.5 Hz
    # =========================
    spikes_05 = r"C:\Internproj\combine\NEURON-C-fiber-main\Results\spikes_Prot22_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"
    stim_05 = r"C:\Internproj\combine\NEURON-C-fiber-main\Results\stim_Prot22_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"

    # =========================
    # 0.6 Hz
    # =========================
    spikes_06 = r"C:\Internproj\combine\NEURON-C-fiber-main\Results\spikes_Prot24_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"
    stim_06 = r"C:\Internproj\combine\NEURON-C-fiber-main\Results\stim_Prot24_gPump-0.0047891_gNav130.010664_gNav170.10664_gNav180.24271_gNav199.48e-05_gKs0.0069733_gKf0.012756_gH0.0025377_gKdr0.018002_gKna0.00042_vRest-55_sineFalse_ampSine0.1.csv"

    aps_04 = pd.read_csv(spikes_04)
    stim_04_df = pd.read_csv(stim_04)

    aps_05 = pd.read_csv(spikes_05)
    stim_05_df = pd.read_csv(stim_05)

    aps_06 = pd.read_csv(spikes_06)
    stim_06_df = pd.read_csv(stim_06)

    plotRecoveryCompareFreq(
        aps_04, stim_04_df,
        aps_05, stim_05_df,
        aps_06, stim_06_df
    )
