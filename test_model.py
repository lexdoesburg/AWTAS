import numpy as np
import matplotlib.pyplot as plt
import time as time_module

import model
import data as data_class

# Test Parameters 1: Page 11 AWTAS
p0 = 3.6e6 # Pa
h = 100 # m
r = 0.05 # m
qm = -0.005 # m^3/s
k = 1e-12 # m^2
phi = 0.1
rho = 813.37 # Water at 240 degrees celsius
nu = 0.0001111 # Water at 240 degrees celsius
C = 0.001303 # Water at 240 degrees celsius
t = 43200 # seconds
time = np.linspace(0, 43200, num=100)
parameters = [p0, qm, h, rho, nu, C, r]

# Build the model and data
theis_model = model.Theis_Solution()
theis_model.generate_data(phi, k, time, parameters, noise=True, sd=1e-3, save_file=True, filename="theis_test1_datafile.txt")
theis_model.data.generate_datafile("theis_test1_datafile2.txt")
opt_phi, opt_k = theis_model.find_model_parameters()
print(opt_phi)
print(opt_k)


# plt.plot(t, data,"k-",label="Synthetic Data (Theis Solution W/Noise)")
start = time_module.time()
plt.semilogx(np.log(time), theis_model.data.observation,"kx",label="Synthetic Data (Theis Solution W/Noise)")
plt.semilogx(np.log(time), theis_model.data.approximation,"r-",label="Approximated Curve")
# plt.plot(t, p, "g-", label="Theis Analytic Solution")
plt.title("Observed Data vs Fitted Curve")
plt.xlabel("Log Time (s)")
plt.ylabel("Pressure (Pa)")
plt.legend(loc="best")
end = time_module.time()
print('Time elapsed = {}'.format(end - start))
plt.show()

# # Test 2 parameters
# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.015 # m^3/s
# k = 4e-12 # m^2
# phi = 0.01
# rho = 527.59 # Water at 360 degrees celsius
# nu = 0.0000603 # Water at 360 degrees celsius
# C = 0.03748 # Water at 360 degrees celsius
# t = 54000 # seconds
# time = np.linspace(0, t, num=100)
# parameters = [p0, qm, h, rho, nu, C, r]

# # Build the model and data
# theis_model = model.Theis_Solution()
# theis_model.generate_data(phi, k, time, parameters, noise=True, sd=2e-3, save_file=True, filename="theis_test2_datafile.txt")
# theis_model.data.generate_datafile("theis_test2_datafile2.txt")
# opt_phi, opt_k = theis_model.find_model_parameters()
# print(opt_phi)
# print(opt_k)

# # plt.plot(t, data,"k-",label="Synthetic Data (Theis Solution W/Noise)")
# start = time_module.time()
# plt.semilogx(np.log(time), theis_model.data.observation,"kx",label="Synthetic Data (Theis Solution W/Noise)")
# plt.semilogx(np.log(time), theis_model.data.approximation,"r-",label="Approximated Curve")
# # plt.plot(t, p, "g-", label="Theis Analytic Solution")
# plt.title("Observed Data vs Fitted Curve")
# plt.xlabel("Log Time (s)")
# plt.ylabel("Pressure (Pa)")
# plt.legend(loc="best")
# end = time_module.time()
# print('Time elapsed = {}'.format(end - start))
# plt.show()


# # Test 3 parameters
# p0 = 3.6e6 # Pa
# h = 100 # m
# r = 0.05 # m
# qm = -0.015 # m^3/s
# k = 1e-14 # m^2
# phi = 0.01
# rho = 527.59 # Water at 360 degrees celsius
# nu = 0.0000603 # Water at 360 degrees celsius
# C = 0.03748 # Water at 360 degrees celsius
# t = 54000 # seconds
# time = np.linspace(0, t, num=100)
# parameters = [p0, qm, h, rho, nu, C, r]

# # Build the model and data
# theis_model = model.Theis_Solution()
# theis_model.generate_data(phi, k, time, parameters, noise=True, sd=1e-8, save_file=True, filename="theis_test3_datafile.txt")
# theis_model.data.generate_datafile("theis_test3_datafile2.txt")
# opt_phi, opt_k = theis_model.find_model_parameters()
# print(opt_phi)
# print(opt_k)

# # plt.plot(t, data,"k-",label="Synthetic Data (Theis Solution W/Noise)")
# start = time_module.time()
# plt.semilogx(np.log(time), theis_model.data.observation,"kx",label="Synthetic Data (Theis Solution W/Noise)")
# plt.semilogx(np.log(time), theis_model.data.approximation,"r-",label="Approximated Curve")
# # plt.plot(t, p, "g-", label="Theis Analytic Solution")
# plt.title("Observed Data vs Fitted Curve")
# plt.xlabel("Log Time (s)")
# plt.ylabel("Pressure (Pa)")
# plt.legend(loc="best")
# end = time_module.time()
# print('Time elapsed = {}'.format(end - start))
# plt.show()