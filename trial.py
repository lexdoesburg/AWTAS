# variables = {}
# values = [{'Value':0.1, 'Units':'xD'}, {'Value':1e-12, 'Units':'D:'}]
# variables['Porosity'], variables['Permeability'] = values
# print(variables)

# print(variables['Porosity']['Value'])
dict_value = dict.fromkeys(['Value', 'Units'])
parameters = {'Initial Pressure' : dict_value, 'Mass Flowrate' : dict_value, 'Layer Thickness' : dict_value, 'Density' : dict_value,
              'Kinematic Viscosity' : dict_value, 'Compressibility' : dict_value, 'Radius' : dict_value}
# print(parameters)

# {'Value' : None, 'Units' : 'Pa'}
# {'Value' : None, 'Units' : 'Kg/s'}
# {'Value' : None, 'Units' : 'm'}
# {'Value' : None, 'Units' : 'Kg/m3'}
# {'Value' : None, 'Units' : 'm2/s'}
# {'Value' : None, 'Units' : '1/Pa'}
# {'Value' : None, 'Units' : 'm'}

# parameters=[None]*7

# if parameters[0]:
#     print(True)
# else:
#     print(False)

for i, (a, b) in enumerate(parameters.items()):
    print(i)
    print(a)
    print(b)
    print(b['Units'])
    # print(a[0])
    # print(a[1]['Units'])
