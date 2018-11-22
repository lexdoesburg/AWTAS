from theis_solution import*
import matplotlib.pyplot as plt

# Page 11 AWTAS
p0 = 3.6e6 # Pa
h = 100 # m
r = 0.05 # m
qm = 0.005 # m^3/s
k = 10e-13 # m^2
phi = 0.1
rho = 813.37 # Water at 240 degrees celsius
nu = 0.0001111 # Water at 240 degrees celsius
C = 0.001303 # Water at 240 degrees celsius
t = np.linspace(1, 43200, num=100)

p = theis_solution(p0, qm, k, h, phi, rho, nu, C, r, t)

# noise_sd = 0.005
# p += p * noise_sd * np.random.randn(p.shape[0])
print(t)
print(p)
