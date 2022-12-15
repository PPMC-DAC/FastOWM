import os
import time

print "Start : %s" % time.ctime()
# Wsize = [2,3,4,5,6,7,8]
ncores = [1,4,8,12,16,20,24]
start = time.time()
for nc in ncores:
    os.system('echo " "')
    os.system('echo " "')
    os.system('echo " "')
    os.system('echo " "')
    time.sleep(20)
    os.system("./cesga_func_OWM.o ../datos/V21_group1_densified_point_cloud 3 20 0.8 %d 5 " % (nc))
print "End : %s" % time.ctime()
end = time.time()
print("Total Execution time: %f hours" % ((end - start)/3600))
