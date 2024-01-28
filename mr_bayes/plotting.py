import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import spines
import os
import re
import sys

dir=sys.argv[1]
iter_file = sys.argv[2]
sampling_freq=float(sys.argv[3])


s = [re.search("^[0-9]*$",i) for i in os.listdir(dir)]
iter_folders = [file for i,file in enumerate(os.listdir(dir)) if s[i]]
iter_folders = sorted(iter_folders)
num_folders= len(iter_folders)

aucs=list()
asdsfs=list()

iterations = np.loadtxt(iter_file)

for i in iter_folders:
    if "asdsfs" in os.listdir(os.path.join(dir,str(i))):
        asdsfs.append(np.loadtxt(os.path.join(dir,str(i),"asdsfs")))
    if "aucs" in os.listdir(os.path.join(dir,str(i))):
        aucs.append(np.loadtxt(os.path.join(dir,str(i),"aucs")))

aucs =  np.transpose(np.array(aucs))
asdsfs = np.transpose(np.array(asdsfs))
print(aucs.shape)
print(asdsfs.shape)

mean_aucs = np.mean(aucs,axis=1)
mean_asdsfs = np.mean(asdsfs,axis=1) 

aucs_25_quantile= np.quantile(aucs,0.25,axis=1)
aucs_75_quantile= np.quantile(aucs,0.75,axis=1)

#asdsfs_25_quantile = np.quantile(asdsfs,0.25,axis=1)
#asdsfs_75_quantile = np.quantile(asdsfs,0.75,axis=1)
asdsfs_25_quantile = np.array([np.quantile(a[a>0],0.25) for a in asdsfs])
asdsfs_75_quantile = np.array([np.quantile(a[a>0],0.75) for a in asdsfs])

matplotlib.style.use('ggplot')
fig, (ax1,ax3) = plt.subplots(1,2)

aucs_lims = [0.5,1.02]
min_auc,max_auc = aucs_lims
max_asdsf = max(asdsfs_75_quantile)
# We want the 1e-2 threshold for ASDSFs to match the 0.8 threshold for our algorithm
# i.e. (max_auc - 0.8) / (max_auc - min_auc) = (log(max_asdsf) - log(1e-2)) / (log(max_asdsf) - log(min_asdsf) )
# or min_asdsf = max_asdsf * ((log(max_asdsf) - log(1e-2)) * (max_auc - min_auc) / (max_auc - 0.8))
min_asdsf = max_asdsf * np.exp( - ((np.log(max_asdsf) - np.log(1e-2)) * (max_auc - min_auc) / (max_auc - 0.8)) )
asdsfs_lims = [min_asdsf,max_asdsf]

color_asdsfs = 'maroon'
ax1.set_xlabel('Iterations')
ax1.set_ylabel('ASDSF', color=color_asdsfs)
ax1.set_ylim(asdsfs_lims)
print(max(iterations)//sampling_freq)
x = np.linspace(sampling_freq,max(iterations),int(max(iterations)//sampling_freq))
ax1.loglog(x, mean_asdsfs, color=color_asdsfs)
ax1.loglog(x, asdsfs_25_quantile,color=color_asdsfs,linestyle = "dashed")
ax1.loglog(x, asdsfs_75_quantile,color=color_asdsfs,linestyle = "dashed")
ax1.plot([sampling_freq,max(iterations)],[1e-2,1e-2],"k--")
ax1.tick_params(axis='y', labelcolor=color_asdsfs)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color_aucs = 'tab:blue'
ax2.set_ylabel('AUC', color=color_aucs)  # we already handled the x-label with ax1
ax2.set_ylim(aucs_lims)
ax2.semilogx(iterations,mean_aucs, color=color_aucs)
ax2.semilogx(iterations, aucs_25_quantile,color=color_aucs,linestyle = "dashed")
ax2.semilogx(iterations, aucs_75_quantile,color=color_aucs,linestyle = "dashed")
ax2.tick_params(axis='y', labelcolor=color_aucs)
ax2.spines["right"].set_color(color_aucs)
ax2.spines["left"].set_color(color_asdsfs)

fig.tight_layout()  

# The iteration file contains multiples of sampling_freq
inds = np.array(list(map(lambda x: int(x//sampling_freq-1),iterations)))

mean_asdsfs = mean_asdsfs[inds]

plt.subplot(1,2,2)
plt.scatter(mean_asdsfs,mean_aucs,c=iterations,norm=matplotlib.colors.LogNorm(),cmap="bone",edgecolor="k")
plt.plot([1e-2,1e-2],aucs_lims,"--",c=color_asdsfs,linewidth=1)
plt.plot(asdsfs_lims,[0.8,0.8],"--",c=color_aucs,linewidth=1)
plt.colorbar()
plt.xscale("log")
plt.ylim(aucs_lims)
plt.xlabel("ASDSF")
# plt.ylabel("AUC")
ax3.spines["bottom"].set_color(color_asdsfs)
ax3.spines["left"].set_color(color_aucs)
ax3.xaxis.label.set_color(color_asdsfs)
ax3.yaxis.label.set_color(color_aucs)
ax3.tick_params(colors=color_asdsfs, axis='x')  
ax3.tick_params(colors=color_aucs, axis='y')  
plt.tight_layout()

plt.show()
