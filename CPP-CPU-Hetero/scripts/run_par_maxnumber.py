import os
import time
import numpy as np

#Sliding window size
Wsize = 10
#Gris size for stage 3 (filling empty areas of Bsize*Bsize with a minimum point)
Bsize = 20
#Overlap for the sliding window. Displacement will be Wsize(1-Overlap)=10m*0.2=2m
Overlap = 0.8
#num_threads for the openmp implementations of stage1 and stage3
num_threads = [1, 8]
#number of times the OWM is executed
nreps = 5
#minRadius for the OWM
mR=0.1

start = time.time()
print("Start : %s" % time.ctime())

executable_par="../bin/parallelmaxnum"
inputs=["../bin/data/Alcoy",
        "../bin/data/Arzua",
        "../bin/data/BrionF",
        "../bin/data/BrionU"]

output="parallelcpp.out"
maxNumber=[32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536]
for file in inputs:
    for mN in maxNumber:
        for nth in num_threads:
            print("Running: {} {} {} {} {} {} {} {} {}".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps,mR,mN))
            f = open(output, "a")
            f.write("\n\nRunning: {} {} {} {} {} {} {} {} {}\n\n".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps,mR,mN))
            f.close()
            os.system("%s %s %d %d %f %d %d %f %d| tee -a %s" % (executable_par, file, Wsize, Bsize, Overlap, nth, nreps,mR, mN, output))

end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
