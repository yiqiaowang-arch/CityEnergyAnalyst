from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
import numpy as np

problem = {
    'num_vars': 5,
    'names': ['Kf', 'Kb', 'b','Kc', 'c'],
    'bounds': [[-0.02178,0.02178],[-0.00099, 0.00099],[0,1],[-0.00099,0.00099],[0,1]]
    #'bounds': [[-2.5,0.0453],[-0.125, 0.125],[0,1],[-0.125,0.125],[0,1]]


}

param_values = saltelli.sample(problem, 35)

print param_values

np.savetxt(r'C:\Users\Zhongming\Desktop\D_test1.csv', param_values, delimiter=',')

