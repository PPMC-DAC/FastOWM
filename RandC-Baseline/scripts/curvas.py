import numpy as np
import matplotlib.pyplot as plt
import sys
import statistics
import os, sys
from matplotlib.backends.backend_pdf import PdfPages
from operator import itemgetter

#Configuration variables
titlefs = 28
ylabelfs = 26
xlabelfs = 20
xticksfs = 20
yticksfs = 28
legendfs = 28
linew = 2
markers = 10

fig1 = plt.figure(figsize=(15, 8))

marks=['o-','x-','s-','+-','v-']
configuration=['1CPU','2CPU','3CPU','4CPU','5CPU','6CPU','7CPU','8CPU']

# plt.plot([2,3,4,5,6,7,8,9,10], [56.71,56.32,57.16,58.12,59.15,60.08,60.76,61.92,62.36], 'o-', linewidth=linew, markersize=markers)
# plt.plot([2,3,4,5,6,7,8,9,10], [55.98,56.88,57.96,58.91,59.75,60.58,61.30,61.99,62.09], 'o-', linewidth=linew, markersize=markers)
plt.plot([2,3,4,5,6,7,8,9,10], [96.23,97.44,98.58,99.52,100.066,100.97,101.51,102.1,103.15], 's-', linewidth=linew, markersize=markers,color='g')
plt.plot([2,3,4,5,6,7,8,9,10], [96,97,98.03,99.31,99.66,100.03,100.59,101.17,101.22], 'o:', linewidth=linew, markersize=12, fillstyle='none', color='b')


sizes=['SVM','NoSVM']

plt.title('Speedup HD 1400 frames',  fontweight='bold', fontsize=titlefs)
plt.legend(sizes,loc='best', fontsize= legendfs)
plt.ylabel('Speedup', fontsize=ylabelfs)
plt.xlabel('Configuration', fontsize=xlabelfs)
#plt.xticks(np.arange(1,len(labels)+1),labels,fontsize=xticksfs)
#Este es para mostrar solo el GPU y el heterogeneo. El anterior es el general para todo.
plt.xticks(np.arange(2,len(configuration)+2),configuration,fontsize=xticksfs)
plt.yticks(fontsize=yticksfs)
plt.grid()

plt.show()
pp = PdfPages("./HD_SVMvsNOSVM.pdf")
pp.savefig(fig1)
pp.close()
