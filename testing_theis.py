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
time = 43200 # seconds
t = np.linspace(0, 43200, num=100)

p = theis_solution(p0, qm, k, h, phi, rho, nu, C, r, t)

# noise_sd = 0.005
# p += p * noise_sd * np.random.randn(p.shape[0])
print(t)
print(p)

# Test if generating data and synthetic data works
data = generate_data(phi, k, 100, time, p0, qm, h, rho, nu, C, r)
noisey_data = generate_data(phi, k, 100, time, p0, qm, h, rho, nu, C, r, noise=True, sd=1e-3)
print(data)
print(noisey_data)
# f, ax1 = plt.subplots(nrows=1, ncols=1)
# ax1.plot(t, data, "k-", label="Theis Data")
# ax1.plot(t, noisey_data, "g-", label="Synthetic Data")
# ax1.set_xlabel("Time (s)")
# ax1.set_ylabel("Pressure (Pa)")
# ax1.legend(loc="best")
# plt.show()

# Test if non-linear optimisation working
# phi = input("Enter porosity: ")
# k = input("Enter permeability: ")
phi, k = find_model_parameters(noisey_data, p0, qm, h, rho, nu, C, r, t, phi=2., k=2.3e-12)
print(phi)
print(k)

# Test how well the parameters fit the data
approximated_data = generate_data(phi, k, 100, time, p0, qm, h, rho, nu, C, r)

plt.plot(t, noisey_data,"k-",label="Synthetic Data (Theis Solution W/Noise)")
plt.plot(t, approximated_data,"r-",label="Approximated Curve")
plt.plot(t, p, "g-", label="Theis Analytic Solution")
plt.title("Observed Data vs Fitted Curve")
plt.xlabel("Time (s)")
plt.ylabel("Pressure (Pa)")
plt.legend(loc="best")

plt.show()