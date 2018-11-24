from model import*
import matplotlib.pyplot as plt
import time as time_module

# Page 11 AWTAS
p0 = 3.6e6 # Pa
h = 100 # m
r = 0.05 # m
qm = -0.005 # m^3/s
k = 1e-12 # m^2
phi = 0.1
rho = 813.37 # Water at 240 degrees celsius
nu = 0.0001111 # Water at 240 degrees celsius
C = 0.001303 # Water at 240 degrees celsius
time = 43200 # seconds
t = np.linspace(0, 43200, num=100)

theis_model = Theis_Solution(p0, qm, h, rho, nu, C, r, t)
theis_model.generate_data(phi, k, noise=True, sd=1e-4, save_file=True)
opt_phi, opt_k = theis_model.find_model_parameters()
actual_data = theis_model.data
print(opt_phi)
print(opt_k)

# Test how well the parameters fit the data
approximated_data = theis_model.generate_data(opt_phi, opt_k)

# plt.plot(t, data,"k-",label="Synthetic Data (Theis Solution W/Noise)")
start = time_module.time()
plt.plot(t, noisey_data,"kx",label="Synthetic Data (Theis Solution W/Noise)")
plt.plot(t, approximated_data,"r-",label="Approximated Curve")
plt.plot(t, p, "g-", label="Theis Analytic Solution")
plt.title("Observed Data vs Fitted Curve")
plt.xlabel("Time (s)")
plt.ylabel("Pressure (Pa)")
plt.legend(loc="best")
end = time_module.time()
print('Time elapsed = {}'.format(end - start))
plt.show()