import os
import time

print "Start : %s" % time.ctime()
ncores = [2,4,6,8,10,12,14,16,18,20,22,24]
start = time.time()
os.system("./serie.o ../datos/V21_group1_densified_point_cloud 10 20 0.8 1 15 ")
for nc in ncores:
    os.system('echo " "')
    os.system('echo " "')
    os.system('echo " "')
    os.system('echo " "')
    time.sleep(20)
    os.system("./parallel.o ../datos/V21_group1_densified_point_cloud 10 20 0.8 %d 15 " % (nc))
print "End : %s" % time.ctime()
end = time.time()
print("Total Execution time: %f hours" % ((end - start)/3600))
