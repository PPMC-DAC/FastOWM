import os
import numpy
import time

nreps = 5

start = time.time()
print("Start : %s" % time.ctime())

sequential="../bin/sequential"
parallel="../bin/parallelgpu"
experiment="gpu"

# for rad in numpy.linspace(0.3,0.4,nreps):
#     output="slurm-Alcoy-"+experiment+"-"+str(rad)+".out"
#     os.system("./ex_alcoy.sh %s %s %f | tee %s" % (sequential,parallel,rad,output))
# for rad in numpy.linspace(0.13,0.2,nreps):
#     output="slurm-Arzua-"+experiment+"-"+str(rad)+".out"
#     os.system("./ex_arzua.sh %s %s %f | tee %s" % (sequential,parallel,rad,output))
# for rad in numpy.linspace(0.1,0.2,nreps):
#     output="slurm-Brion-forestal-"+experiment+"-"+str(rad)+".out"
#     os.system("./ex_brionf.sh %s %s %f | tee %s" % (sequential,parallel,rad,output))
for rad in numpy.linspace(0.1,0.2,nreps):
    output="slurm-Brion-urban-"+experiment+"-"+str(rad)+".out"
    os.system("./ex_brion.sh %s %s %f | tee %s" % (sequential,parallel,rad,output))


# sequential="../bin/sequential"
# parallel="../bin/task"
# experiment="medSize-tr"

# for rad in numpy.linspace(0.9,1.2,nreps):
#     output="slurm-Alcoy-"+experiment+"-"+str(rad)+".out"
#     os.system("./ex_alcoy.sh %s %s %f | tee %s" % (sequential,parallel,rad,output))
# for rad in numpy.linspace(0.8,1.1,nreps):
#     output="slurm-Arzua-"+experiment+"-"+str(rad)+".out"
#     os.system("./ex_arzua.sh %s %s %f | tee %s" % (sequential,parallel,rad,output))
# for rad in numpy.linspace(0.2,0.5,nreps):
#     output="slurm-Brion-forestal-"+experiment+"-"+str(rad)+".out"
#     os.system("./ex_brionf.sh %s %s %f | tee %s" % (sequential,parallel,rad,output))
# for rad in numpy.linspace(0.2,0.5,nreps):
#     output="slurm-Brion-urban-"+experiment+"-"+str(rad)+".out"
#     os.system("./ex_brion.sh %s %s %f | tee %s" % (sequential,parallel,rad,output))

end = time.time()
print("End : %s" % time.ctime())
print("Total Execution time: %f hours" % ((end - start)/3600))