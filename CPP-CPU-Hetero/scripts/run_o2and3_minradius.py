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

start = time.time()
print("Start : %s" % time.ctime())

executable_par="../bin/o2and3_minrad"
inputs=["../bin/data/AlcoyH",
        "../bin/data/ArzuaH",
        "../bin/data/BrionFH",
        "../bin/data/BrionUH"]

output="o2and3_minradius.out"
minRadius=list(np.arange(0.1,2,0.1))

for file in inputs:
    for mR in minRadius:
        for nth in num_threads:
            print("Running: {} {} {} {} {} {} {} {}".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps,mR))
            f = open(output, "a")
            f.write("\n\nRunning: {} {} {} {} {} {} {} {}\n\n".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps,mR))
            f.close()
            os.system("%s %s %d %d %f %d %d %f| tee -a %s" % (executable_par, file, Wsize, Bsize, Overlap, nth, nreps,mR, output))

end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
