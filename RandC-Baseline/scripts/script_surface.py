
from matplotlib.backends.backend_pdf import PdfPages
from operator import itemgetter
from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt
import sys
import statistics
import os, sys

# Fichero de entrada
input_text=sys.argv[1]
output_text=sys.argv[2]

os.system("txt2las -i %s -o %s --parse xyz" % (input_text,output_text))


# def read_values(input_text):
#     values=[]
#     with open(input_text) as f:
#         for line in f:
#             # Divido la linea en tantos elementos como tenga
#             tokens = line.split()
#             values.append((tokens))
#             # if len(tokens) and tokens[0] in keys[0]:
#             #     if(int(tokens[2])==sizes[0]):
#             #         values=get_time_values(tokens,sizes[0],f)
#     return values




# print("SCRIPT PYTHON 3")
# values = read_values(input_text)
# values=np.array(values)
# values = values.astype(np.float)
# print("Número de puntos a pintar ",len(values))
# # print(values[:,2])
#
# X, Y, Z = values[:,0], values[:,1], values[:,2]
#
# fig = plt.figure()
# ax = plt.axes(projection='3d')
#
# # Data for three-dimensional scattered points
# zdata = Z
# xdata = Y
# ydata = X
# # ax.plot_trisurf(xdata, ydata, zdata, color='white', edgecolors='grey', alpha=0.2)
# ax.scatter3D(xdata, ydata, zdata, color='green');
#
# plt.show()
