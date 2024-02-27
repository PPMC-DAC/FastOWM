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
#num_threads for the openmp implementations of stage1 and stage3
num_threads = 1
#number of times the OWM is executed
nreps = 1

# name of the input files without the extension .xyz
inputs=[
    "../bin/data/AlcoyH",
    "../bin/data/ArzuaH",
    "../bin/data/BrionFH",
    "../bin/data/BrionUH",
    ]

# get the hostname
hostname = os.popen("hostname").read().strip()
# executables
executable="../bin/sequential"

start = time.time()
print("Start : %s" % time.ctime())

for cloud in inputs:
    # set the output file
    output = f'{cloud}-seq_{hostname}.out'
    print("Running: {} {} {} {} {} {} {}".format(executable,cloud,Wsize,Bsize,Overlap,1,nreps))
    os.system("%s %s %d %d %f %d %d | tee %s" % (executable, cloud, Wsize, Bsize, Overlap, 1, nreps, output))


end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
