import numpy as np
import matplotlib.pyplot as plt
import time as time_module

import model_copy as model
import data as data_class


"""
Determined that leastsq is most effective and is equivalent to an unbounded curve_fit
"""


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
t = 43200 # seconds
time = np.linspace(0, 43200, num=100)
parameters = [p0, qm, h, rho, nu, C, r]

# Print actual values
print('Actual phi: {} k: {}\n'.format(phi, k))

# Build the model and data


# # Test leastsq
# theis_model1 = model.Theis_Solution()
# theis_model1.generate_data(phi, k, time, parameters, noise=True, sd=1e-3)
# opt_phi, opt_k = theis_model1.find_model_parameters()
# print('leastsq phi: {} k: {}\n'.format(opt_phi, opt_k))

# # Test curve_fit
# theis_model2 = model.Theis_Solution()
# theis_model2.generate_data(phi, k, time, parameters, noise=True, sd=1e-3)
# opt_phi, opt_k = theis_model2.find_model_parameters(curvefit=True)
# print('curvefit phi: {} k: {}\n'.format(opt_phi, opt_k))

# # Test least_squares
# theis_model3 = model.Theis_Solution()
# theis_model3.generate_data(phi, k, time, parameters, noise=True, sd=1e-3)
# opt_phi, opt_k = theis_model3.find_model_parameters(leastsquares=True)
# print('least_squares phi: {} k: {}\n'.format(opt_phi, opt_k))

# 46773860381.72343
# 652074936.5173881


# # plt.plot(t, data,"k-",label="Synthetic Data (Theis Solution W/Noise)")
# start = time_module.time()
# plt.semilogx(np.log(time), theis_model1.data.observation,"kx",label="Synthetic Data (Theis Solution W/Noise)")
# plt.semilogx(np.log(time), theis_model1.data.approximation,"r-",label="leastsq")
# plt.semilogx(np.log(time), theis_model2.data.approximation,"g-",label="curve_fit")
# plt.semilogx(np.log(time), theis_model3.data.approximation,"b--",label="least_squares")
# # plt.plot(t, p, "g-", label="Theis Analytic Solution")
# plt.title("Observed Data vs Fitted Curve")
# plt.xlabel("Log Time (s)")
# plt.ylabel("Pressure (Pa)")
# plt.legend(loc="best")
# end = time_module.time()
# print('Time elapsed = {}'.format(end - start))
# plt.show()

estimated_phi = []
estimated_k = []
optimal_phi_leastsq = []
optimal_k_leastsq = []
optimal_phi_curvefit = []
optimal_k_curvefit = []
optimal_phi_least_squares = []
optimal_k_least_squares = []
optimal_phi_curvefit_bound = []
optimal_k_curvefit_bound = []
optimal_phi_least_squares_bound = []
optimal_k_least_squares_bound = []

ls_error = []
cf_error = []
cfb_error = []
lsq_error = []
lsqb_error = []

theis_model1 = model.Theis_Solution()
theis_model1.generate_data(phi, k, time, parameters, noise=True, sd=1e-3)

for phi1 in np.linspace(0.01, 0.2, 20):
    for k1 in [1e-16, 0.5e-16, 1e-15, 0.5e-15, 1e-14, 0.5e-14, 1e-13, 0.5e-13, 1e-12, 0.5e-12]:
        estimated_phi.append(phi1)
        estimated_k.append(k1)

        # Test leastsq
        opt_phi, opt_k = theis_model1.find_model_parameters(phi=phi1, k=k1)
        optimal_phi_leastsq.append(opt_phi)
        optimal_k_leastsq.append(opt_k)
        ls_error.append(theis_model1.weighted_error())
        
        # Test curve_fit unbounded
        opt_phi, opt_k = theis_model1.find_model_parameters(phi=phi1, k=k1, curvefit=True)
        optimal_phi_curvefit.append(opt_phi)
        optimal_k_curvefit.append(opt_k)
        cf_error.append(theis_model1.weighted_error())

        # Test least_squares unbounded 
        opt_phi, opt_k = theis_model1.find_model_parameters(phi=phi1, k=k1, leastsquares=True)
        optimal_phi_least_squares.append(opt_phi)
        optimal_k_least_squares.append(opt_k)
        cfb_error.append(theis_model1.weighted_error())

        # Test curve_fit bounded
        opt_phi, opt_k = theis_model1.find_model_parameters(phi=phi1, k=k1, curvefit=True, bounded=True)
        optimal_phi_curvefit_bound.append(opt_phi)
        optimal_k_curvefit_bound.append(opt_k)
        lsq_error.append(theis_model1.weighted_error())

        # Test least_squares bounded 
        opt_phi, opt_k = theis_model1.find_model_parameters(phi=phi1, k=k1, leastsquares=True, bounded=True)
        optimal_phi_least_squares_bound.append(opt_phi)
        optimal_k_least_squares_bound.append(opt_k)
        lsqb_error.append(theis_model1.weighted_error())

# with open('nonlinear_solver.csv', 'w') as file:
#     file.write('phi guess, k guess, phi leastsq, k leastsq, phi curve_fit, k curve_fit, phi least_squares, k least_squares, phi curve_fit bounded, k curve_fit bounded, phi least_squares bounded, k least_squares bounded\n')
#     for i in range(len(estimated_phi)):
#         file.write('{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(estimated_phi[i], estimated_k[i], optimal_phi_leastsq[i], optimal_k_leastsq[i], optimal_phi_curvefit[i], optimal_k_curvefit[i], optimal_phi_least_squares[i], optimal_k_least_squares[i], optimal_phi_curvefit_bound[i], optimal_k_curvefit_bound[i], optimal_phi_least_squares_bound[i], optimal_k_least_squares_bound[i]))

with open('nl_error.txt', 'w') as file:
    file.write('leastsq, curve_fit, curve_fit bounded, least_squares, least_squares bounded\n')
    for i in range(len(ls_error)):
        file.write('{}, {}, {}, {}, {}\n'.format(ls_error[i], cf_error[i], cfb_error[i], lsq_error[i], lsqb_error[i]))

    