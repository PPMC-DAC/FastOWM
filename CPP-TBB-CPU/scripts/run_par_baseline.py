import os
import time
import numpy as np

# Alder
nprocs = 16 # 8 p-cores and 8 e-cores

#Sliding window size
Wsize = 10
#Gris size for stage 3 (filling empty areas of Bsize*Bsize with a minimum point)
Bsize = 20
#Overlap for the sliding window. Displacement will be Wsize(1-Overlap)=10m*0.2=2m
Overlap = 0.8
#num_threads for the openmp implementations of stage1 and stage3
# num_threads = [1, 2, 4, 6, 8]
num_threads = np.insert(np.linspace(2, nprocs, nprocs//2, dtype=int, endpoint=True), 0, 1)
print(num_threads)
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

for file in inputs:
    output="baseline.out"
    print("Running: {} {} {} {} {} {} {}".format(executable_seq,file,Wsize,Bsize,Overlap,1,nreps))
    f = open(output, "a")
    f.write("\n\nRunning: {} {} {} {} {} {} {}\n\n".format(executable_seq,file,Wsize,Bsize,Overlap,1,nreps))
    f.close()
    os.system("%s %s %d %d %f %d %d | tee -a %s" % (executable_seq, file, Wsize, Bsize, Overlap, 1, nreps, output))
    for nth in num_threads:
        print("Running: {} {} {} {} {} {} {}".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps))
        f = open(output, "a")
        f.write("\n\nRunning: {} {} {} {} {} {} {}\n\n".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps))
        f.close()
        os.system("%s %s %d %d %f %d %d | tee -a %s" % (executable_par, file, Wsize, Bsize, Overlap, nth, nreps, output))

end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
