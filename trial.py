import numpy as np

data = data_class.Data()
new_model = model.SKG9D(data)
new_model.generate_data(variables, time, parameters=None, noise = True, sd = 300, save_file=True, filename='SKG9D_test1.dat')