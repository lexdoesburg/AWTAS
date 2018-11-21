from theis_solution import*
import matplotlib.pyplot as plt

# Page 11 AWTAS
p0 = 3.6e6 # Pa
h = 100 # m
r = 0.5 # m
q = 0.005 # m^3/s
k = 10e-13 # m^2
porosity = 0.1
density = 813.37 # Water at 240 degrees celsius
v = 0.0001111 # Water at 240 degrees celsius
C = 0.001303 # Water at 240 degrees celsius
t = np.linspace(1, 43200, num=100)

p = theis_solution(p0, q, k, h, porosity, density, v, C, r, t)

noise_sd = 0.005
p += p * noise_sd * np.random.randn(p.shape[0])

f, ax1 = plt.subplots(nrows=1, ncols=1)

ax1.plot(t/3600, p/1e5)
ax1.set_xlabel("Time (hours)")
ax1.set_ylabel("Pressure (bar)")
# testing github
#again
plt.show()