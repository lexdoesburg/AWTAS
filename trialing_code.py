from theis_solution import*
import matplotlib.pyplot as plt

# x = np.linspace(0,4)
# e1 = exp1(x)
# f, ax1 = plt.subplots(nrows=1, ncols=1)

# ax1.plot(x,e1)
# ax1.set_xlabel("x")
# ax1.set_ylabel("e1")

# plt.show()

phi=0.1
k=10e-13

x0 = phi, k
print(x0)

initial_parameters = np.array([phi, k])
a, b = initial_parameters
print(a)
print(b)