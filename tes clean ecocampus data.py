import pandas as pd
import matplotlib.pylab as plt
import numpy as np
from os import listdir
from os.path import isfile, join

mypath = r'C:\Users\JimenoF\Documents\GitHub\cea-reference-case\reference-case-ecocampus\baseline\inputs\building-metering/'
path_demand = r'C:\Users\JimenoF\Documents\GitHub\cea-reference-case\reference-case-ecocampus\baseline\outputs\data\demand/Total_demand.csv'
# buildings = [f for f in listdir(mypath) if isfile(join(mypath, f))]
#B081
cluster = 'C004.csv'
buildings_in_cluster = ['B201',
'B202'
]
areas = pd.read_csv(path_demand, usecols=['Name', 'GFA_m2']).set_index('Name')


n = len(buildings_in_cluster)


new_areas = areas.T[buildings_in_cluster]
total_area = new_areas.ix['GFA_m2'].sum()


for building in buildings_in_cluster:
    data = pd.read_csv(mypath + '//' + cluster, usecols=['Name', 'Ef_kWh'])
    wighterd_average =  new_areas.ix['GFA_m2', building]/ total_area
    data['Ef_kWh'] = data['Ef_kWh']*wighterd_average
    data['Name'] = building
    data.to_csv(mypath+building+'.csv')


# for name in buildings:
#     fileraw = pd.read_csv(mypath+ name)
#
#     time =  fileraw['DATE']
#     tale =  fileraw['Ef_kWh'][:8752].values
#     beginning = fileraw['Ef_kWh'][8752:8760].values
#     new_values= np.append(beginning, tale)
#
#     # fixing range
#     fileraw['Ef_kWh'] = new_values
#     fileraw.to_csv(mypath+ name)
#     # fileraw.plot()
#     # plt.show()
