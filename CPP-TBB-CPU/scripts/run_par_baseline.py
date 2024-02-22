import os
import time
import numpy as np

# get the number of physical cores
nprocs = int(os.popen("lscpu | grep 'Core(s) per socket' | awk '{print $4}'").read().strip())

#Sliding window size
Wsize = 10
#Gris size for stage 3 (filling empty areas of Bsize*Bsize with a minimum point)
Bsize = 20
#Overlap for the sliding window. Displacement will be Wsize(1-Overlap)=10m*0.2=2m
Overlap = 0.8
#num_threads for the openmp implementations of stage1 and stage3
num_threads = np.insert(np.linspace(2, nprocs, nprocs//2, dtype=int, endpoint=True), 0, 1)
#number of times the OWM is executed
nreps = 5

start = time.time()
print("Start : %s" % time.ctime())

executable_seq="../bin/baseline"
executable_par="../bin/par_baseline"
inputs=["../bin/data/AlcoyH",
        "../bin/data/ArzuaH",
        "../bin/data/BrionFH",
        "../bin/data/BrionUH"]

# get the hostname
hostname = os.popen("hostname").read().strip()
# set the output file
output = f'baseline_{hostname}.out'

with open(output, "a") as f:
    for cloud in inputs:
        print("Running: {} {} {} {} {} {} {}".format(executable_seq,cloud,Wsize,Bsize,Overlap,1,nreps))
        # save the configuration in the file
        f.write("\n\nRunning: {} {} {} {} {} {} {}\n\n".format(executable_seq,cloud,Wsize,Bsize,Overlap,1,nreps))
        # flush the buffer
        f.flush()
        # execute the command and save the output to the file
        os.system("%s %s %d %d %f %d %d | tee -a %s" % (executable_seq, cloud, Wsize, Bsize, Overlap, 1, nreps, output))
        for nth in num_threads:
            print("Running: {} {} {} {} {} {} {}".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps))
            # save the configuration in the file
            f.write("\n\nRunning: {} {} {} {} {} {} {}\n\n".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps))
            # flush the buffer
            f.flush()
            # execute the command and save the output to the file
            os.system("%s %s %d %d %f %d %d | tee -a %s" % (executable_par, cloud, Wsize, Bsize, Overlap, nth, nreps, output))

end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
# copy the output file to the results folder
os.system(f'cp {output} ../../Results/baseline_{hostname}.txt')
