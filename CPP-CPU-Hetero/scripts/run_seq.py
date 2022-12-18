import os
import numpy
import time

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

start = time.time()
print("Start : %s" % time.ctime())

executable="../bin/sequential"
inputs=["../bin/data/Alcoy",
        "../bin/data/Arzua",
        "../bin/data/BrionF",
        "../bin/data/BrionU"]


for file in inputs:
    output=file+"-seq.out"
    print("Running: {} {} {} {} {} {} {}".format(executable,file,Wsize,Bsize,Overlap,num_threads,nreps))
    os.system("%s %s %d %d %f %d %d | tee %s" % (executable, file, Wsize, Bsize, Overlap, num_threads, nreps, output))



end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))