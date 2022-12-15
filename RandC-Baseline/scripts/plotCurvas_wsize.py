import numpy as np
import matplotlib.pyplot as plt
import sys
import statistics
import os, sys
from matplotlib.backends.backend_pdf import PdfPages
from operator import itemgetter

keys = []
keys.append('BEST:')

values = []
size = [4,8,16,20]
# size = [2,3,4,5,6,7,8]
# size = [1,2,3,4,5,6,7,8]
# size = [0.65,0.7,0.75,0.8,0.85,0.9]

keys.append('BEST:')

def read_values(input_text):
    values = []
    with open(input_text) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) and tokens[0] in keys[0]:
                # print("Lectura de BEST: %f" % float(tokens[1]))
                values.append(float(tokens[1]))

    return values

input="V21_inter"
input2="V21_inter2"
input3="VaihingenS09_MAX6"
input4="VaihingenS09_MAX8"
input5="VaihingenS09_MAX10"
# input="slurm_Guitiriz.out"
# fig1 = plt.figure(1)
values=read_values(input+".out")
values=np.array(values)
# input="slurm_Guitiriz_foto.out"
# fig2 = plt.figure(2)
values1=read_values(input2+".out")
values1=np.array(values1)
# input="slurm_Begonte.out"
# fig3 = plt.figure(3)
values2=read_values(input3+".out")
values2=np.array(values2)

values3=read_values(input4+".out")
values3=np.array(values3)

values4=read_values(input5+".out")
values4=np.array(values4)

# # input="slurm_Guitiriz.out"
# fig1 = plt.figure(1)
# # values=read_values(input)
# # values=np.array(values)
# plt.bar(cores,values, color="lightblue")
# plt.plot(cores,values, 'o-', linewidth=2, markersize=6, color="green")
# plt.title('TIME ON FLY: Guitiriz',  fontweight='bold', fontsize=12)
# plt.ylabel('Time', fontsize=12)
# plt.xlabel('Cores', fontsize=12)
# plt.xticks(cores,fontsize=10)
# plt.yticks(fontsize=10)
# plt.ylim(0,5)
# # print("Se han encontrado %d experimentos\n" % len(values))
#
# # input="slurm_Guitiriz_foto.out"
# fig2 = plt.figure(2)
# # values1=read_values(input)
# # values1=np.array(values1)
# plt.bar(cores,values1, color="lightblue")
# plt.plot(cores,values1, 'o-', linewidth=2, markersize=6, color="green")
# plt.title('TIME ON FLY: Guitiriz_foto',  fontweight='bold', fontsize=12)
# plt.ylabel('Time', fontsize=12)
# plt.xlabel('Cores', fontsize=12)
# plt.xticks(cores,fontsize=10)
# plt.yticks(fontsize=10)
# plt.ylim(0,9)
#
# # input="slurm_Begonte.out"
# fig3 = plt.figure(3)
# # values2=read_values(input)
# # values2=np.array(values2)
# plt.bar(cores,values2, color="lightblue")
# plt.plot(cores,values2, 'o-', linewidth=2, markersize=6, color="green")
# plt.title('TIME ON FLY: Begonte',  fontweight='bold', fontsize=12)
# plt.ylabel('Time', fontsize=12)
# plt.xlabel('Cores', fontsize=12)
# plt.xticks(cores,fontsize=10)
# plt.yticks(fontsize=10)
# plt.ylim(0,14)

# fig7 = plt.figure(7)
# plt.plot(cores,np.log10(values), 'o-', linewidth=2, markersize=6, color="green")
# plt.plot(cores,np.log10(values1), 'o-', linewidth=2, markersize=6, color="blue")
# plt.plot(cores,np.log10(values2), 'o-', linewidth=2, markersize=6, color="orange")
# plt.title('TIME ON FLY',  fontweight='bold', fontsize=12)
# plt.ylabel('Time', fontsize=12)
# plt.xlabel('Cores', fontsize=12)
# plt.xticks(cores,fontsize=10)
# plt.yticks(fontsize=10)
# plt.legend(["Guitiriz","Guitiriz_foto","Begonte"],loc='best', fontsize=16)
# # plt.ylim(0,14)
# plt.grid()

fig4 = plt.figure(4)
plt.plot(size,values, 'o-', linewidth=2, markersize=6, color="green")
plt.plot([2,4,8,16,20,24],values1, 'o-', linewidth=2, markersize=6, color="blue")
# plt.plot(size,values2, 'o-', linewidth=2, markersize=6, color="orange")
# plt.plot(size,values3, 'o-', linewidth=2, markersize=6)
# plt.plot(size,values4, 'o-', linewidth=2, markersize=6)
plt.title('TIME ON FLY',  fontweight='bold', fontsize=12)
plt.ylabel('Time [s]', fontsize=12)
plt.xlabel('Cores', fontsize=12)
plt.xticks(size,fontsize=10)
plt.yticks(fontsize=10)
plt.legend(["W=2 B=2 O=0.9","W=12 B=20 O=0.9"],loc='best', fontsize=16)
# plt.ylim(0,14)
plt.grid()

# speedup=values[0]/values
# fig5 = plt.figure(5)
# plt.plot(cores,speedup,'o-', linewidth=2, markersize=6, color="green")
# speedup=values1[0]/values1
# plt.plot(cores,speedup,'o-', linewidth=2, markersize=6, color="blue")
# speedup=values2[0]/values2
# plt.plot(cores,speedup,'o-', linewidth=2, markersize=6, color="orange")
# plt.title('SPEEDUP',  fontweight='bold', fontsize=12)
# plt.ylabel('Improvement [x]', fontsize=12)
# plt.xlabel('Cores', fontsize=12)
# plt.xticks(cores,fontsize=10)
# plt.yticks(fontsize=10)
# plt.legend([input,input2,input3],loc='best', fontsize=16)
# plt.grid()

# efficiency=values[0]/(cores*values)
# fig6 = plt.figure(6)
# plt.plot(cores,efficiency,'o-', linewidth=2, markersize=6, color="green")
# efficiency=values1[0]/(cores*values1)
# plt.plot(cores,efficiency,'o-', linewidth=2, markersize=6, color="blue")
# efficiency=values2[0]/(cores*values2)
# plt.plot(cores,efficiency,'o-', linewidth=2, markersize=6, color="orange")
# plt.title('EFFICIENCY',  fontweight='bold', fontsize=12)
# plt.ylabel('Performance', fontsize=12)
# plt.xlabel('Cores', fontsize=12)
# plt.xticks(cores,fontsize=10)
# plt.yticks(fontsize=10)
# plt.legend([input,input2,input3],loc='best', fontsize=16)
# plt.grid()

# speedup=values[0]/values
# fig5 = plt.figure(5)
# plt.plot(cores,speedup,'o-', linewidth=2, markersize=6, color="green")
# z = np.polyfit(cores,speedup,0.5)
# p = np.poly1d(z)
# print("len fit: %d" %len(p(cores)))
# plt.plot(cores,np.array(p(cores)),'o:', linewidth=2, markersize=6, color="green")
# plt.title('SPEEDUP: Guitiriz',  fontweight='bold', fontsize=12)
# plt.ylabel('Improvement', fontsize=12)
# plt.xlabel('Cores', fontsize=12)
# plt.xticks(cores,fontsize=10)
# plt.yticks(fontsize=10)
# plt.grid()

plt.show()
