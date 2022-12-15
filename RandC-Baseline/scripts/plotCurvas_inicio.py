import numpy as np
import matplotlib.pyplot as plt
import sys
import statistics
import os, sys
from matplotlib.backends.backend_pdf import PdfPages
from operator import itemgetter

keys= []
values=[]
num_execution_per_input=5
# sizes = [500,600,700,800]
sizes = [7000]
configuration=['GPU','1CPU\nGPU','2CPU\nGPU','3CPU\nGPU','4CPU\nGPU','5CPU\nGPU','6CPU\nGPU','7CPU\nGPU','8CPU\nGPU']
# input_text=os.path.splitext(sys.argv[1])[0]

keys.append('BEST:')
# keys.append('Time:')
# keys.append('Frames')



# Path to be created
# path = "./"+input_text

# os.makedirs(path, exist_ok=True);

#print "Path is created"




def read_values(input_text):
    values = []
    with open(input_text) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) and tokens[0] in keys[0]:
                print("Lectura de BEST: %f" % float(tokens[1]))
                values.append(float(tokens[1]))
                # if(int(tokens[2])==sizes[0]):
                #     values=get_time_values(tokens,sizes[0],f)
                # elif(int(tokens[2])==sizes[1]):
                #     values=get_time_values(tokens,sizes[1],f)
                # elif(int(tokens[2])==sizes[2]):
                #     values=get_time_values(tokens,sizes[2],f)
                # elif(int(tokens[2])==sizes[3]):
                #    values=get_time_values(tokens,sizes[3],f)
                #elif(int(tokens[2])==sizes[4]):
                #    values=get_time_values(tokens,sizes[4],f)

    # values.sort(key=itemgetter(1,0))
    return values

def get_time_values(tokens,size,f):
    next_line1 = f.readline()
    tokens = next_line1.split()
    token_value=float(tokens[1])
    myfile = open("./"+path+"/Speedup_Values_"+sys.argv[1], 'w')
    if(int(tokens[4])== 0):
        next_line = f.readline()
        token_frame = next_line.split()
        cpu_value=float(token_frame[6])
        gpu_value=float(token_frame[2])
        total_frames=float(token_frame[8])
        if(total_frames!=size):
            print("Error found: Not computed all frames in execution "+str(next_line1))
            print("Error found: Not computed all frames in execution "+str(next_line1),file=myfile)
        values.append((2,size,token_value,gpu_value,cpu_value))
    if(int(tokens[4])== 1):
        if(int(tokens[7])==0):
            next_line = f.readline()
            token_frame = next_line.split()
            cpu_value=float(token_frame[6])
            gpu_value=float(token_frame[2])
            total_frames=float(token_frame[8])
            if(total_frames!=size):
                print("Error found: Not computed all frames in execution"+str(next_line1))
                print("Error found: Not computed all frames in execution "+str(next_line1),file=myfile)
            values.append((1,size,token_value,gpu_value,cpu_value))
        else:
            next_line = f.readline()
            token_frame = next_line.split()
            cpu_value=float(token_frame[6])
            gpu_value=float(token_frame[2])
            total_frames=float(token_frame[8])
            if(total_frames!=size):
                print("Error found: Not computed all frames in execution"+str(next_line1))
                print("Error found: Not computed all frames in execution "+str(next_line1),file=myfile)
            values.append((3,size,token_value,gpu_value,cpu_value))
    if(int(tokens[4])== 2):
        next_line = f.readline()
        token_frame = next_line.split()
        cpu_value=float(token_frame[6])
        gpu_value=float(token_frame[2])
        total_frames=float(token_frame[8])
        if(total_frames!=size):
            print("Error found: Not computed all frames in execution"+str(next_line1))
            print("Error found: Not computed all frames in execution "+str(next_line1),file=myfile)
        values.append((4,size,token_value,gpu_value,cpu_value))
    if(int(tokens[4])== 3):
        next_line = f.readline()
        token_frame = next_line.split()
        cpu_value=float(token_frame[6])
        gpu_value=float(token_frame[2])
        total_frames=float(token_frame[8])
        if(total_frames!=size):
            print("Error found: Not computed all frames in execution"+str(next_line1))
            print("Error found: Not computed all frames in execution "+str(next_line1),file=myfile)
        values.append((5,size,token_value,gpu_value,cpu_value))
    if(int(tokens[4])== 4):
        next_line = f.readline()
        token_frame = next_line.split()
        cpu_value=float(token_frame[6])
        gpu_value=float(token_frame[2])
        total_frames=float(token_frame[8])
        if(total_frames!=size):
            print("Error found: Not computed all frames in execution"+str(next_line1))
            print("Error found: Not computed all frames in execution "+str(next_line1),file=myfile)
        values.append((6,size,token_value,gpu_value,cpu_value))
    if(int(tokens[4])== 5):
        next_line = f.readline()
        token_frame = next_line.split()
        cpu_value=float(token_frame[6])
        gpu_value=float(token_frame[2])
        total_frames=float(token_frame[8])
        if(total_frames!=size):
            print("Error found: Not computed all frames in execution"+str(next_line1))
            print("Error found: Not computed all frames in execution "+str(next_line1),file=myfile)
        values.append((7,size,token_value,gpu_value,cpu_value))
    if(int(tokens[4])== 6):
        next_line = f.readline()
        token_frame = next_line.split()
        cpu_value=float(token_frame[6])
        gpu_value=float(token_frame[2])
        total_frames=float(token_frame[8])
        if(total_frames!=size):
            print("Error found: Not computed all frames in execution"+str(next_line1))
            print("Error found: Not computed all frames in execution "+str(next_line1),file=myfile)
        values.append((8,size,token_value,gpu_value,cpu_value))
    if(int(tokens[4])== 7):
        next_line = f.readline()
        token_frame = next_line.split()
        cpu_value=float(token_frame[6])
        gpu_value=float(token_frame[2])
        total_frames=float(token_frame[8])
        if(total_frames!=size):
            print("Error found: Not computed all frames in execution"+str(next_line1))
            print("Error found: Not computed all frames in execution "+str(next_line1),file=myfile)
        values.append((9,size,token_value,gpu_value,cpu_value))
    if(int(tokens[4])== 8):
        next_line = f.readline()
        token_frame = next_line.split()
        cpu_value=float(token_frame[6])
        gpu_value=float(token_frame[2])
        total_frames=float(token_frame[8])
        if(total_frames!=size):
            print("Error found: Not computed all frames in execution"+str(next_line1))
            print("Error found: Not computed all frames in execution "+str(next_line1),file=myfile)
        values.append((10,size,token_value,gpu_value,cpu_value))
    return(values)

def get_min_time(values):
    times=[]
    #values.remove(min(values))
    #for i in range(int(len(values)/2)+1):
    #    values.remove(max(values))
    #times.append(statistics.mean(values))
    times.append(float(min(values)))
    return times

def get_final_times(values):
    y=[]
    x=[]
    z=[]
    time=[]
    #time.clear()
    count=0
    #for i in range(len(values)):
    for j in range(1,len(values)+1):
        if(j%num_execution_per_input!=0):
            y.append(values[j-1][2])
            x.append(values[j-1][3])
            z.append(values[j-1][4])
            #count+=1
            #print(j)
        else:
            y.append(values[j-1][2])
            x.append(values[j-1][3])
            z.append(values[j-1][4])
            #print(y)
            time.append((get_min_time(y),values[j-1][0],values[j-1][1],get_min_time(x),get_min_time(z)))
            #print(values[j])
            y.clear()
            x.clear()
            z.clear()
            #count=0
    return time



def plot_speedup(times,sizes,labels):
    sp = []
    ef =[]
    values=[]
    threads=[]
    for i in range(len(times)):
        threads.append(times[i][1])
        values.append(times[i][0])
    values=np.array(values)
    threads=np.array(threads)
    j=0
    sp = []
    myfile = open("./"+path+"/Speedup_Values_"+sys.argv[1], 'w')
    print(threads-2)
    #print(values)
    for (v,x) in zip(values,threads):
        if int(x) == 1:
            j=j+1
            ser=float(v)
        print("Cores: %s \t Sec: %s \t Speed-up: %s" % (x-2, v, ser/float(v)))
        print("Cores: %s \t Sec: %s \t Speed-up: %s" % (x-2, v, ser/float(v)),file=myfile)
        spd=ser/float(v)
        sp.append(spd)
        ef.append(spd/int(x))


    puntos=int(len(sp)/j)
    print ("Número de experimentos a pintar: %d " % j)
    print ("Número de puntos en cada experimento: %d" % puntos)
    print ("Número de experimentos a pintar: %d " % j,file=myfile)
    print ("Número de puntos en cada experimento: %d" % puntos,file=myfile)
    myfile.close()

    #Configuration variables
    titlefs = 20
    ylabelfs = 18
    xlabelfs = 18
    xticksfs = 16
    yticksfs = 16
    legendfs = 14
    linew = 2
    markers = 8

    fig1 = plt.figure(figsize=(15, 8))

    marks=['o-','x-','s-','+-','v-']
    configuration=['1GPU','1CPU\nGPU','2CPU\nGPU','3CPU\nGPU','4CPU\nGPU','5CPU\nGPU','6CPU\nGPU','7CPU\nGPU','8CPU\nGPU']
    for i in range(j):
        init=(i)*puntos+1# el +4 es para saltarse la parte multicore y empezar con el 1GPU en el speedup
        end= (i+1)*(puntos)
       # print np.array(threads[init:end]), np.array(sp[init:end])
        plt.plot(np.array(threads[init:end]), np.array(sp[init:end]), marks[i], linewidth=linew, markersize=markers)

    #plt.plot(np.array(threads[4:9]), np.array(sp[4:9]), marks[1], linewidth=linew, markersize=markers)
    #plt.plot(np.array(threads[13:18]), np.array(sp[13:18]), marks[2], linewidth=linew, markersize=markers)
    #plt.plot(np.array(threads[22:27]), np.array(sp[22:27]), marks[3], linewidth=linew, markersize=markers)
    #plt.plot(np.array(threads[0:puntos]), np.array(threads[0:puntos]), '-', linewidth=linew, markersize=markers)

    #sizes.append('Ideal')

    plt.title('Speedup',  fontweight='bold', fontsize=titlefs)
    plt.legend(sizes,loc='best', fontsize= legendfs)
    plt.ylabel('Speedup', fontsize=ylabelfs)
    plt.xlabel('Configuration', fontsize=xlabelfs)
    #plt.xticks(np.arange(1,len(labels)+1),labels,fontsize=xticksfs)
    #Este es para mostrar solo el GPU y el heterogeneo. El anterior es el general para todo.
    plt.xticks(np.arange(2,len(configuration)+2),configuration,fontsize=xticksfs)
    plt.yticks(fontsize=yticksfs)
    plt.grid()

    plt.show()
    pp = PdfPages("./"+path+"/Speedup-"+sys.argv[1]+".pdf")
    pp.savefig(fig1)
    pp.close()
'''
    configuration=['1CPU-GPU','2CPU-GPU','3CPU-GPU','4CPU-GPU','5CPU-GPU','6CPU-GPU']

    tiempo_cpu=[]
    tiempo_gpu=[]
    for i in range(len(times)):
        tiempo_cpu.append(times[i][4])
        tiempo_gpu.append(times[i][3])
    tiempo_cpu=np.array(tiempo_cpu)
    tiempo_gpu=np.array(tiempo_gpu)

    for j in range(len(sizes)):
        fig2 = plt.figure(figsize=(15,8))
        fig2.suptitle('Workload distribution CPU-GPU for size :'+str(sizes[j]),fontweight='bold', fontsize=16)
        for i in range(6):
            cpu=float(tiempo_cpu[i+2+j*9])
            gpu=float(tiempo_gpu[i+2+j*9])
            # cpu=float(times[i+2][4])
            # gpu=float(times[i+2][3])
            # print ("Porcentaje CPU: %s " % times[i+2][4])
            # print ("Porcentaje GPU: %s " % times[i+2][3])
            # Data to plot
            labels = 'CPU', 'GPU'
            # devices = [times[i+1+(j*5)][4], times[i+1+(j*5)][3]]
            devices=[cpu,gpu]
            colors = ['lightcoral', 'lightskyblue']
            explode = (0.1, 0)  # explode 1st slice
            plt.subplot(2, 3,i+1)
            plt.title(configuration[i], fontsize=titlefs,fontweight='bold')
            # Plot
            plt.pie(devices,explode=explode, labels=labels, colors=colors,autopct='%1.1f%%', startangle=100,textprops={'fontsize': 18})#,'fontweight':'bold'
            plt.legend(labels, loc="best",fontsize=titlefs)
            plt.axis('equal')

        plt.show()
        pp = PdfPages("./"+path+"/CPU-GPU_"+str(sizes[j])+"_"+input_text+".pdf")
        pp.savefig(fig2)
        pp.close()

'''

cores=range(1,25)
input="slurm_Guitiriz.out"
fig1 = plt.figure(1)
values=read_values(input)
plt.bar(cores,values, color="lightblue")
plt.plot(cores,values, 'o-', linewidth=2, markersize=6, color="green")
print("Se han encontrado %d experimentos\n" % len(values))


input="slurm_Begonte.out"
fig2 = plt.figure(2)
values=read_values(input)
plt.bar(cores,values, color="lightblue")
plt.plot(cores,values, 'o-', linewidth=2, markersize=6, color="green")
print("Se han encontrado %d experimentos\n" % len(values))
plt.show()
