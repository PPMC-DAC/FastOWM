import os
import time
from datetime import datetime
from common.utils import get_nprocs


# name of the input files without the extension .xyz
inputs=[
        "../bin/data/INAER_2011_Alcoy_Core.xyz",
        "../bin/data/BABCOCK_2017_Arzua_3B.xyz",
        "../bin/data/V21_group1_densified_point_cloud.xyz",
        "../bin/data/V19_group1_densified_point_cloud.xyz",
        ]

# get the hostname
hostname = os.popen("hostname").read().strip()
# set the output file
output = f"sycl_cuda_{hostname}.out"
#num_threads for the openmp implementations of stage1 and stage3
num_threads = get_nprocs()

# executables
executable_list = [ 
                    "../bin/owm-sycl-cpu",
                    "../bin/owm-sycl-cpu-nomemo",
                    "../bin/owm-sycl-igpu",
                    "../bin/owm-sycl-igpu-nomemo",
                    "../bin/owm-sycl-dgpu",
                    "../bin/owm-sycl-dgpu-nomemo",
                    "../bin/owm-cuda",
                    "../bin/owm-cuda-grid",
                    "../bin/owm-cuda-nomemo",
                ]

maxNumber=[4,8,16,32,64,128,256,512,1024]

# select the number of threads
if hostname == 'coffeelake1':
    vnth = [8,16] # physical and hyperthreading
elif hostname == 'alder':
    # only p-cores, and p-cores + e-cores
    vnth = [8,12,16,24]
elif hostname == 'bombay':
    vnth = [32,48,64,80,96]
else:
    # select the maximum number of threads
    vnth = [num_threads[-1]]

start = time.time()
print("Start : %s" % time.ctime())

with open(output, "a") as f:
    # stamp the start time
    f.write(f'Start: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}')
    # iterate over the clouds, levels and number of threads
    for exe in executable_list:
        for i in inputs:
            for mN in maxNumber:
                if 'cpu' in exe:
                    # if cpu version, iterate over the number of threads
                    for nth in vnth:
                        print("\n***************\nRunning: {} {} {} {}".format(exe, i, mN, nth))
                        # save the configuration in the file
                        f.write("\n\nRunning: {} {} {} {}\n\n".format(exe, i, mN, nth))
                        # flush the buffer
                        f.flush()
                        # execute the command and save the output to the file
                        os.system("DPCPP_CPU_NUM_CUS=%d %s %s %d | tee -a %s" % (nth, exe, i, mN, output))
                        # sleep until the next execution
                        time.sleep(30)
                else:
                    print("\n***************\nRunning: {} {} {}".format(exe, i, mN))
                    # save the configuration in the file
                    f.write("\n\nRunning: {} {} {}\n\n".format(exe, i, mN))
                    # flush the buffer
                    f.flush()
                    # execute the command and save the output to the file
                    os.system("%s %s %d | tee -a %s" % (exe, i, mN, output))
                    # sleep until the next execution
                    time.sleep(30)

    end = time.time()
    f.write(f'End: {datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}\n')
    f.write("Total Execution time: {} hours\n".format((end - start)/3600))
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
# copy the output file to the results folder
os.system(f"cp {output} ../../Results/{output.replace('.out', '.txt')}")


