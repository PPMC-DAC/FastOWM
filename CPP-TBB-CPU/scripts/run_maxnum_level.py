import os
import time
import numpy as np
from datetime import datetime
from common.utils import get_nprocs

#Sliding window size
Wsize = 10
#Gris size for stage 3 (filling empty areas of Bsize*Bsize with a minimum point)
Bsize = 20
#Overlap for the sliding window. Displacement will be Wsize(1-Overlap)=10m*0.2=2m
Overlap = 0.8
#num_threads for the openmp implementations of stage1 and stage3
num_threads = get_nprocs()
#number of times the OWM is executed
nreps = 5

# name of the input files without the extension .xyz
inputs=[
    "../bin/data/AlcoyH",
    "../bin/data/ArzuaH",
    "../bin/data/BrionFH",
    "../bin/data/BrionUH",
    ]

# get the hostname
hostname = os.popen("hostname").read().strip()
# set the output file
output = f'o4_maxnumber_{hostname}.out'
# executables
executable_par="../bin/o4maxnum"
# results file in .csv format
resultsFile = output.replace('.out', '.csv')

maxNumber=[8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536]
levels = list(range(3,10))

# select the number of threads
if hostname == 'alder':
    # only p-cores, and p-cores + e-cores
    vnth = [8,12,16,24]
elif hostname == 'bombay':
    vnth = [32,48,64,80,96]
else:
    # select the maximum number of threads
    vnth = [num_threads[-1]]

mR=0.1

start = time.time()
print("Start : %s" % time.ctime())

with open(output, "a") as f:
    # stamp the start time
    f.write(f'Start: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}')
    # iterate over the clouds, levels and number of threads
    for cloud in inputs:
        for mN in maxNumber:
            for lev in levels:
                for nth in vnth:
                    print("Running: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps,mR,mN,lev))
                    # save the configuration in the file
                    f.write("\n\nRunning: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}\n\n".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps,mR,mN,lev))
                    # flush the buffer
                    f.flush()
                    # excecute the command and save the output to the file
                    os.system("%s -i %s -W %d -B %d -O %f -n %d -l %d -r %f -s %d -L %d --results %s| tee -a %s" % (executable_par, cloud, Wsize, Bsize, Overlap, nth, nreps, mR, mN, lev, resultsFile, output))

    end = time.time()
    f.write(f'End: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}\n')
    f.write("Total Execution time: {} hours\n".format((end - start)/3600))

print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
# copy the output file to the results folder
os.system(f"cp {output} ../../Results/{output.replace('.out', '.txt')}")
os.system(f"cp {resultsFile} ../../Results/{resultsFile}")
