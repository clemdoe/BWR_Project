### IMPORTS

import sys
import matplotlib.pyplot as plt

### MAIN

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

ax1.set_ylabel("power (W)")
ax2.set_ylabel("reactivity (pcm)")
ax3.set_ylabel("fuelTemperature (K)")
ax3.set_xlabel("time(s)")

for q in range(len(sys.argv)-1) :

    filename = sys.argv[q+1]
    file = open(filename, "r")
    lines = file.readlines()
    powers = []
    times = []
    totRhos = []
    TFuels = []
    i = 0
    power = 0
    time = 0
    totRho = 0
    TFuel = 0
    for line in lines :
        if "totalPower" in line :
            power = float(line.split()[2])
        if "totalReactivity" in line :
            totRho = float(line.split()[2])
        if "TFuel = " in line :
            TFuel = float(line.split()[2])
        if "Time =" in line and not "ExecutionTime" in line :
            time = float(line.split()[2])
            powers.append(power)
            times.append(time)
            totRhos.append(totRho)
            TFuels.append(TFuel)
    del totRhos[0]
    del powers[0]
    del TFuels[0]
    t0 = times[0]
    for i in range(len(times)) :
    	times[i] = times[i] - t0  
    del times[-1]
    ax1.semilogx(times, powers, label=filename)
    ax2.semilogx(times, totRhos, label=filename)  
    ax3.semilogx(times, TFuels, label=filename)      
ax1.legend()
ax1.grid(True)
ax2.legend()
ax2.grid(True)
ax3.legend()
ax3.grid(True)
plt.show()