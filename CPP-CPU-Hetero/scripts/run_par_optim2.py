import os
import time

#Sliding window size
Wsize = 10
#Gris size for stage 3 (filling empty areas of Bsize*Bsize with a minimum point)
Bsize = 20
#Overlap for the sliding window. Displacement will be Wsize(1-Overlap)=10m*0.2=2m
Overlap = 0.8
#number of times the OWM is executed
nreps = 5

output="o2parallel.out"
start = time.time()
print("Start : %s" % time.ctime())
f = open(output, "a")
f.write("Start : {}".format(time.ctime()))
f.close()

executable_par="../bin/o2par"
inputs=["../bin/data/AlcoyH",
        "../bin/data/ArzuaH",
        "../bin/data/BrionFH",
        "../bin/data/BrionUH"]

num_threads = [1, 2, 4, 6, 8]
levels = list(range(3,9))
mR=0.1

for file in inputs:
    for lev in levels:
        for nth in num_threads:
            print("Running: {} {} {} {} {} {} {} {} {}".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps,mR,lev))
            f = open(output, "a")
            f.write("\n\nRunning: {} {} {} {} {} {} {} {} {}\n\n".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps,mR,lev))
            f.close()
            os.system("%s %s %d %d %f %d %d %f %d | tee -a %s" % (executable_par, file, Wsize, Bsize, Overlap, nth, nreps, mR, lev, output))

end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
f = open(output, "a")
f.write("End : {}".format(time.ctime()))
f.write("Total Execution time: {} hours".format((end - start)/3600))
f.close()