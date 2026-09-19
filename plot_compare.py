# plot comparison for membrane potential and sodium current

import pandas as pd
import matplotlib.pyplot as plt


combined_file = (
    r"C:\Internproj\combine\NEURON-C-fiber-main\Results"
    r"\combined_overall.csv"
)

original_file = (
    r"C:\Internproj\original\NEURON-C-fiber-main\Results"
    r"\original_overall.csv"
)

# Load datas
combined = pd.read_csv(combined_file)
original = pd.read_csv(original_file)


# 1.Vm
###############################
plt.figure()
plt.plot(
    combined["Time"],
    combined["Vm"],
    label="3-channel modified",
    color="tab:pink"
)
plt.plot(
    original["Time"],
    original["Vm"],
    label="Original",
    color="tab:blue"
)


plt.xlabel("time (ms)")
plt.ylabel("Vm (mV)")
plt.title("Membrane potential")

plt.legend()
plt.grid()


# 2.Nav current
###############################
plt.figure()
plt.plot(
    combined["Time"],
    combined["Ina"],
    label="3-channel modified",
    color="tab:pink"
)
plt.plot(
    original["Time"],
    original["Ina"],
    label="Original",
    color="tab:blue"
)


plt.xlabel("time (ms)")
plt.ylabel("current (mA/cm2)")
plt.title("Total sodium current")

plt.legend()
plt.grid()


plt.show()
