import os
import time

#Sliding window size
Wsize = 10
#Gris size for stage 3 (filling empty areas of Bsize*Bsize with a minimum point)
Bsize = 20
#Overlap for the sliding window. Displacement will be Wsize(1-Overlap)=10m*0.2=2m
Overlap = 0.8
#num_threads for the openmp implementations of stage1 and stage3
num_threads = [1]
#number of times the OWM is executed
nreps = 5
maxNumber=65536

output="sycl_statsminrad.out"


executable_par="../bin/stats"
inputs=["../bin/data/AlcoyH",
        "../bin/data/ArzuaH",
        "../bin/data/BrionFH",
        "../bin/data/BrionUH"]

#minRadius=[0.9,0.6,0.2,0.2] #Best minRadius according to tree+owm time without memo (see o2and3_minradius.ipynb)
levels=[5,7,9]
nth=8
mR=0.1

#for file,mR in zip(inputs,minRadius):
for lev in levels:
    start = time.time()
    print("Start : %s" % time.ctime())
    f = open(output, "a")
    f.write("Start : {}".format(time.ctime()))
    f.close()
    for file in inputs:
        print("Running: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps,mR,maxNumber,lev))
        f = open(output, "a")
        f.write("\n\nRunning: {} -i {} -W {} -B {} -O {} -n {} -l {} -r {} -s {} -L {}\n\n".format(executable_par,file,Wsize,Bsize,Overlap,nth,nreps,mR,maxNumber,lev))
        f.close()
        os.system("%s -i %s -W %d -B %d -O %f -n %d -l %d -r %f -s %d -L %d| tee -a %s" % (executable_par, file, Wsize, Bsize, Overlap, nth, nreps,mR, maxNumber, lev, output))

    end = time.time()
    print("End : %s" % time.ctime())
    print("Total Execution time: %f hours" % ((end - start)/3600))
    f = open(output, "a")
    f.write("End : {}".format(time.ctime()))
    f.write("Total Execution time: {} hours".format((end - start)/3600))
    f.close()

    os.system("mkdir ../bin/data/MinRadHisto"+str(lev)+" && mv ../bin/data/*.csv ../bin/data/MinRadHisto"+str(lev)+" && mv sycl_statsminrad.out ../bin/data/MinRadHisto"+str(lev)+"/sycl_statsminrad.txt")

