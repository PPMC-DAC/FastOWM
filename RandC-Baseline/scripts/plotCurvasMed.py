import numpy as np
import matplotlib.pyplot as plt
import sys
import statistics
import os, sys
from matplotlib.backends.backend_pdf import PdfPages
from operator import itemgetter

keys = []
keys.append('BEST:')

values = [0,0,0,0,0,0,0,0,0]
values1 = [0,0,0,0,0,0,0,0,0]
values2 = [0,0,0,0,0,0,0,0,0]
# cores=[1,4,8,12,16,20,24]
cores = [1,2,4,6,8,10,12,14,16]

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

input="slurm_V21py"
# values = read_values(input+".out")
values=np.array(values)
for iter in range(1,7):
    values = values + np.array(read_values(input+str(iter)+".out"))
    print(input+str(iter)+".out")
values=values/9

input="slurm_BABpy"
# values = read_values(input+".out")
values1=np.array(values1)
for iter in range(1,7):
    values1 = values1 + np.array(read_values(input+str(iter)+".out"))
    print(input+str(iter)+".out")
values1=values1/9

input="slurm_Corepy"
# values = read_values(input+".out")
values2=np.array(values2)
for iter in range(1,7):
    values2 = values2 + np.array(read_values(input+str(iter)+".out"))
    print(input+str(iter)+".out")
values2=values2/9

# input="slurm_7AM_V21"
# input2="slurm_7AM_V211"
# input3="slurm_7AM_V212"
# input4="slurm_7AM_V213"
# input="slurm_Guitiriz.out"
# fig1 = plt.figure(1)
# values=read_values(input+".out")
# values=np.array(values)
# input="slurm_Guitiriz_foto.out"
# fig2 = plt.figure(2)
# values1=read_values(input2+".out")
# values1=np.array(values1)
# input="slurm_Begonte.out"
# fig3 = plt.figure(3)
# values2=read_values(input3+".out")
# # values2=np.array([values2[0],values2[2],values2[4],values2[6],values2[8],values2[10],values2[12]])
# values2=np.array(values2)
#
# values3=read_values(input4+".out")
# values3=np.array(values3)

# values=(values+values1+values2+values3)/4

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
plt.plot(cores,values, 'o-', linewidth=2, markersize=6, color="green")
plt.plot(cores,values1, 'o-', linewidth=2, markersize=6, color="blue")
plt.plot(cores,values2, 'o-', linewidth=2, markersize=6, color="orange")
plt.title('TIME ON FLY',  fontweight='bold', fontsize=12)
plt.ylabel('Time [s]', fontsize=12)
plt.xlabel('Cores', fontsize=12)
plt.xticks(cores,fontsize=10)
plt.yticks(fontsize=10)
plt.legend(["Fotogrametría","BABCOCK","Alcoy"],loc='best', fontsize=16)
# plt.ylim(0,14)
plt.grid()

speedup=values[0]/values
fig5 = plt.figure(5)
plt.plot(cores,cores, 'o--', linewidth=2, markersize=6, color="red") # ideal
plt.plot(cores,speedup,'o-', linewidth=2, markersize=6, color="green")
speedup=values1[0]/values1
plt.plot(cores,speedup,'o-', linewidth=2, markersize=6, color="deepskyblue")
speedup=values2[0]/values2
plt.plot(cores,speedup,'o-', linewidth=2, markersize=6, color="purple")
plt.title('SPEEDUP',  fontweight='bold', fontsize=12)
plt.ylabel('Improvement [x]', fontsize=12)
plt.xlabel('Cores', fontsize=12)
plt.xticks(cores,fontsize=10)
plt.yticks(fontsize=10)
plt.legend(["IDEAL","Fotogrametría","BABCOCK","Alcoy"],loc='best', fontsize=16)
plt.grid()

fig6 = plt.figure(6)
plt.plot(cores,np.ones(9), 'o--', linewidth=2, markersize=6, color="red") # ideal
efficiency=values[0]/(cores*values)
plt.plot(cores,efficiency,'o-', linewidth=2, markersize=6, color="green")
efficiency=values1[0]/(cores*values1)
plt.plot(cores,efficiency,'o-', linewidth=2, markersize=6, color="deepskyblue")
efficiency=values2[0]/(cores*values2)
plt.plot(cores,efficiency,'o-', linewidth=2, markersize=6, color="purple")
plt.title('EFFICIENCY',  fontweight='bold', fontsize=12)
plt.ylabel('Performance', fontsize=12)
plt.xlabel('Cores', fontsize=12)
plt.xticks(cores,fontsize=10)
plt.yticks(fontsize=10)
plt.legend(["IDEAL","Fotogrametría","BABCOCK","Alcoy"],loc='best', fontsize=16)
plt.grid()

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
