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
num_threads = [1, 2, 4, 6, 8]
#number of times the OWM is executed
nreps = 5
maxNumber=[32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536]

output="sycl_maxnumber.out"

start = time.time()
print("Start : %s" % time.ctime())
f = open(output, "a")
f.write("Start : {}".format(time.ctime()))
f.close()

executable_par="../bin/parallelmaxnum"
inputs=["../bin/data/AlcoyH",
        "../bin/data/ArzuaH",
        "../bin/data/BrionFH",
        "../bin/data/BrionUH"]

mR=0.1
lev=5

for file in inputs:
    for mN in maxNumber:
        for nth in num_threads:
            print("Running: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps,mR,mN,lev))
            f = open(output, "a")
            f.write("\n\nRunning: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}\n\n".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps,mR,mN,lev))
            f.close()
            os.system("%s -i %s -W %d -B %d -O %f -n %d -l %d -r %f -s %d -L %d| tee -a %s" % (executable_par, file, Wsize, Bsize, Overlap, nth, nreps,mR, mN, lev, output))


end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
f = open(output, "a")
f.write("End : {}".format(time.ctime()))
f.write("Total Execution time: {} hours".format((end - start)/3600))
f.close()