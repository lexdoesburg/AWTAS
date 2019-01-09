import data as datastructure
import numpy as np
import matplotlib.pyplot as plt
import time


# ---------------- Tough2 data
# logtime, barpressure = np.genfromtxt('Pwell.dat', delimiter='   ', skip_header=0).T
# logtime2, barpressure2 = np.genfromtxt('Pwell2.dat', delimiter='   ', skip_header=0).T


# for i in range(len(logtime)):
#     print(logtime[i], barpressure[i])

# plt.plot(logtime,barpressure,'k-')
# plt.plot(logtime2,barpressure2,'r-')

# # plt.semilogx(logtime,barpressure,'k-')
# # plt.semilogx(logtime2,barpressure2,'r-')

# plt.show()
# -----------------
datapoints = []
pressure = []
for i in range(10000):
    p = np.random.rand()
    datapoints.append(datastructure.DataPoint(i, p))
    pressure.append(p)
    print(datapoints[i].time, datapoints[i].observation)

start = time.clock()
pressure2 = [p.observation for p in datapoints]
end = time.clock()
print(end-start)

start = time.clock()
pressure2 = pressure
end = time.clock()
print(end-start)