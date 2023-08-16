import os
import time

output="sycl_output.out"

start = time.time()
print("Start : %s" % time.ctime())
f = open(output, "a")
f.write("Start : {}".format(time.ctime()))
f.close()

executable=["../bin/owm-sycl-cpu","../bin/owm-sycl-cpu-nomemo","../bin/owm-sycl-igpu","../bin/owm-sycl-igpu-nomemo"]
inputs=["../../LiDARClouds/INAER_2011_Alcoy_Core.xyz",
        "../../LiDARClouds/BABCOCK_2017_Arzua_3B.xyz",
        "../../LiDARClouds/V21_group1_densified_point_cloud.xyz",
        "../../LiDARClouds/V19_group1_densified_point_cloud.xyz"]

maxNumber=[4,8,16,32,64,128,256,512,1024]

for exe in executable:
    for i in inputs:
        for mN in maxNumber:
            print("\n***************\nRunning: {} {} {}".format(exe, i, mN))
            f = open(output, "a")
            f.write("\n\nRunning: {} {} {}\n\n".format(exe, i, mN))
            f.close()
            os.system("%s %s %d | tee -a %s" % (exe, i, mN, output))
            time.sleep(30)

end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))
f = open(output, "a")
f.write("End : {}".format(time.ctime()))
f.write("Total Execution time: {} hours".format((end - start)/3600))
f.close()



