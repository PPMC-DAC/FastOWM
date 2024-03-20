import os
import time
import numpy as np
from datetime import datetime
from common.utils import get_nprocs

#Sliding window size
Wsize = 10
#Gris size for stage 3 (filling empty areas of Bsize*Bsize with a minimum point)
Bsize = 20
#Overlap for the sliding window. Displacement will be Wsize(1-Overlap)=10m*0.2=2m
Overlap = 0.8
#num_threads for the openmp implementations of stage1 and stage3
num_threads = get_nprocs()
#number of times the OWM is executed
nreps = 5

# name of the input files without the extension .xyz
inputs=[
    "../bin/data/AlcoyH",
    "../bin/data/ArzuaH",
    "../bin/data/BrionFH",
    "../bin/data/BrionUH",
    ]

# get the hostname
hostname = os.popen("hostname").read().strip()
# set the output file
output_list = [ f'baseline-rev1_{hostname}.out', f'baseline-rev1_collapse_{hostname}.out', # dynamic
                # f'baseline-rev1_static_{hostname}.out', f'baseline-rev1_collapse_static_{hostname}.out', # static
                # f'baseline-rev1_guided_{hostname}.out', f'baseline-rev1_collapse_guided_{hostname}.out', # guided
                # f'baseline-rev1_tasks_{hostname}.out', # tasks
                ]
# executables
executable_list = [ "../bin/baseline-rev1", "../bin/baseline-rev1collap",
                    # "../bin/baseline-rev1static", "../bin/baseline-rev1collapstatic",
                    # "../bin/baseline-rev1guided", "../bin/baseline-rev1collapguided",
                    # "../bin/baseline-rev1tasks"
                    ]

# list of chunk sizes used in the dynamic scheduling of stage 1
# chunk_list = list(range(1,9))
chunk_list = [1, 2, 4, 8, 10, 12, 14]

# zip the output and executable lists
for output, executable_par in zip(output_list, executable_list):

    start = time.time()
    print("Start : %s" % time.ctime())

    with open(output, "a") as f:
        # stamp the start time
        f.write(f'Start: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}')
        # iterate over the clouds, levels and number of threads
        for cloud in inputs:
            for chunk in chunk_list:
                for nth in num_threads:
                    print("Running: {} {} {} {} {} {} {} {}".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps, chunk))
                    # save the configuration in the file
                    f.write("\n\nRunning: {} {} {} {} {} {} {} {}\n\n".format(executable_par,cloud,Wsize,Bsize,Overlap,nth,nreps, chunk))
                    # flush the buffer
                    f.flush()
                    # execute the command and save the output to the file
                    os.system("%s %s %d %d %f %d %d %d | tee -a %s" % (executable_par, cloud, Wsize, Bsize, Overlap, nth, nreps, chunk, output))

        end = time.time()
        f.write(f'End: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}\n')
        f.write("Total Execution time: {} hours\n".format((end - start)/3600))

    print("End : %s" % time.ctime())
    print("Total Execution time: %f hours" % ((end - start)/3600))
    # copy the output file to the results folder
    os.system(f"cp {output} ../../Results/{output.replace('.out', '.txt')}")