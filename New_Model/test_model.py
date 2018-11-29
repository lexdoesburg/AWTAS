import numpy as np
import matplotlib.pyplot as plt
import time as time_module

import model
import data as data_class

# Test New Model
data = data_class.Data()
data.read_file('SKG9D_press.csv')
new_model = model.SKG9D()

print(data.time)
print(data.observation)