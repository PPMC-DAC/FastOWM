import os
import time

print "Start : %s" % time.ctime()
# Wsize = [2,3,4,5,6,7,8]
ncores = [2,4,6,8,10,12,14,16,18,20,22,24]
start = time.time()
os.system("./serie.o ../datos/INAER_2011_Alcoy_Core 10 20 0.8 1 15 ")
for nc in ncores:
    os.system('echo " "')
    os.system('echo " "')
    os.system('echo " "')
    os.system('echo " "')
    time.sleep(20)
    os.system("./parallel.o ../datos/INAER_2011_Alcoy_Core 10 20 0.8 %d 15 " % (nc))
print "End : %s" % time.ctime()
end = time.time()
print("Total Execution time: %f hours" % ((end - start)/3600))
