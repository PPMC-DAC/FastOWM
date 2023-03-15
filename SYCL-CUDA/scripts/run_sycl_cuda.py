import os
import time

output="sycl_cuda_output.out"

start = time.time()
print("Start : %s" % time.ctime())
f = open(output, "a")
f.write("Start : {}".format(time.ctime()))
f.close()

executable=["../bin/owm-sycl-cpu","../bin/owm-sycl-cpu-nomemo","../bin/owm-sycl-igpu","../bin/owm-sycl-igpu-nomemo",
            "../bin/owm-sycl-dgpu","../bin/owm-sycl-dgpu-nomemo","../bin/owm-cuda","../bin/owm-cuda-nomemo"]
inputs=["../bin/data/INAER_2011_Alcoy_Core.xyz",
        "../bin/data/BABCOCK_2017_Arzua_3B.xyz",
        "../bin/data/V21_group1_densified_point_cloud.xyz",
        "../bin/data/V19_group1_densified_point_cloud.xyz"]

maxNumber=[8,16,32,64,128,256,512,1024]

for exe in executable:
    for i in inputs:
        for mN in maxNumber:
            print("\n***************\nRunning: {} {} {}".format(exe, i, mN))
            f = open(output, "a")
            f.write("\n\nRunning: {} {} {}\n\n".format(exe, i, mN))
            f.close()
            os.system("%s %s %d | tee -a %s" % (exe, i, mN, output))

end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
f = open(output, "a")
f.write("End : {}".format(time.ctime()))
f.write("Total Execution time: {} hours".format((end - start)/3600))
f.close()



