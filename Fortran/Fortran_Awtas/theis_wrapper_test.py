from theis_wrapper import theis

P0 = 3.6e6 # Pa
b = 100 # m
r = 0.05 # m
Q0 = -0.005 # m^3/s
k = 1e-12 # m^2
phi = 0.1
rho = 813.37 # Water at 240 degrees celsius
nu = 0.0001111 # Water at 240 degrees celsius
c = 0.001303 # Water at 240 degrees celsius
t0 = 0
dt = 200
t1 = 54000
numData = 271

pressure = theis(k, nu, phi, rho, c, b, Q0, P0, r, t0, dt, t1, numData)
print(pressure)

