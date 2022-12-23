import os
import time

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

start = time.time()
print("Start : %s" % time.ctime())

executable_par="../bin/o1qtree"
inputs=["../bin/data/Alcoy",
        "../bin/data/Arzua",
        "../bin/data/BrionF",
        "../bin/data/BrionU"]

for file in inputs:
    output="qtree.out"
    for nth in num_threads:
        print("Running: {} {} {} {} {} {} {}".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps))
        f = open(output, "a")
        f.write("\n\nRunning: {} {} {} {} {} {} {}\n\n".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps))
        f.close()
        os.system("%s %s %d %d %f %d %d | tee -a %s" % (executable_par, file, Wsize, Bsize, Overlap, nth, nreps, output))

end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
