import os
import time
import numpy as np
from datetime import datetime

#Sliding window size
Wsize = 10
#Gris size for stage 3 (filling empty areas of Bsize*Bsize with a minimum point)
Bsize = 20
#Overlap for the sliding window. Displacement will be Wsize(1-Overlap)=10m*0.2=2m
Overlap = 0.8
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
output = f'sycl_statsminrad_{hostname}.out'
# executables
executable_par="../bin/stats"

mN=65536
# minRadius=[0.9,0.6,0.1,0.1] #Best minRadius according to tree+owm time without memo (see o2and3_minradius.ipynb)
# let's get different minRadius values for each cloud
minRadius = [2.0, 0.9, 0.6, 0.1]
lev=5
nth=1

start = time.time()
print("Start : %s" % time.ctime())

with open(output, "a") as f:
    # stamp the start time
    f.write(f'Start: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}')
    # iterate over the clouds and minRadius
    for cloud in inputs:
        for mR in minRadius:
            print("Running: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps,mR,mN,lev))
            # save the configuration in the file
            f.write("\n\nRunning: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}\n\n".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps,mR,mN,lev))
            # flush the buffer
            f.flush()
            # excecute the command and save the output to the file
            os.system("%s -i %s -W %d -B %d -O %f -n %d -l %d -r %f -s %d -L %d --results %s | tee -a %s" % (executable_par, cloud, Wsize, Bsize, Overlap, nth, nreps, mR, mN, lev, f'{cloud}_hist_mR{mR}.csv', output))

    end = time.time()
    f.write(f'End: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}\n')
    f.write("Total Execution time: {} hours\n".format((end - start)/3600))
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))

results_dir = f'../bin/data/MinRadHisto-{hostname}'
# ensure the directory exists
os.makedirs(results_dir, exist_ok=True)
os.system(f"mv ../bin/data/*.csv {results_dir} && mv {output} {results_dir}/{output.replace('.out', '.txt')}")

# set the output file
output = f'sycl_statsmaxnum_{hostname}.out'
# executable
executable_par="../bin/statsmaxnum"

maxNumber=[1024,512,256,128,64,32] #Best maxNumber according to tree+owm time without memo (see o2and3_maxnumber.ipynb)
mR=0.1
lev=5
nth=1

start = time.time()
print("Start : %s" % time.ctime())

with open(output, "a") as f:
    # stamp the start time
    f.write(f'Start: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}')
    # iterate over the clouds and maxNumber
    for cloud in inputs:
        for mN in maxNumber:
            print("Running: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps,mR,mN,lev))
            # save the configuration in the file
            f.write("\n\nRunning: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}\n\n".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps,mR,mN,lev))
            # flush the buffer
            f.flush()
            # excecute the command and save the output to the file
            os.system("%s -i %s -W %d -B %d -O %f -n %d -l %d -r %f -s %d -L %d --results %s | tee -a %s" % (executable_par, cloud, Wsize, Bsize, Overlap, nth, nreps, mR, mN, lev, f'{cloud}_hist_mN{mN}.csv', output))

    end = time.time()
    f.write(f'End: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}\n')
    f.write("Total Execution time: {} hours\n".format((end - start)/3600))
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))

results_dir = f'../bin/data/MaxNumHisto-{hostname}'
# ensure the directory exists
os.makedirs(results_dir, exist_ok=True)
os.system(f"mv ../bin/data/*.csv {results_dir} && mv {output} {results_dir}/{output.replace('.out', '.txt')}")
