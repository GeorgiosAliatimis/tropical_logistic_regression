import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import re

#dirs = [f"phyml{i}" for i in range(10)]

#dirs = [f"primates_4_{i}" for i in [4000,2000,1000]] 
dirs=["primates_1"]
#dirs = [f"primates_{i}_1000" for i in [4,10]] 

# dirs = ["primates_10_1000","primates_4_1000"]
# dirs = ["primates_15000_4_1000","primates_10_1000"]

# dirs = ["six_primates_1000"]
# dirs = ["primates_15000_4_1000"]

#dirs.append("primates_10_60_4000")

#dirs =  ["five_primates_10_1000"]
#dirs = ["phyml0_5000"]

fig = plt.figure()
ax = fig.add_subplot(111)

for dir in dirs:
    s = [re.search("^[0-9]*$",i) for i in os.listdir(dir)]
    iter_folders = [file for i,file in enumerate(os.listdir(dir)) if s[i]]
    iter_folders = sorted(iter_folders)
    num_folders= len(iter_folders)

    num_folders = 100
    iter_folders = iter_folders[:num_folders]

    # aucs=[None]* num_folders
    aucs_tlr=[]
    asdsfs=[]

    for ind,i in enumerate(iter_folders):
        if "asdsfs" in os.listdir(os.path.join(dir,str(i))):
            asdsfs.append(np.loadtxt(os.path.join(dir,str(i),"asdsfs")))
        if "aucs_0.3" in os.listdir(os.path.join(dir,str(i))):
            aucs_tlr.append(np.loadtxt(os.path.join(dir,str(i),"aucs_0.3")))

    # aucs = np.transpose(np.array(aucs))
    aucs_tlr =  np.transpose(np.array(aucs_tlr))
    asdsfs = np.transpose(np.array(asdsfs))
    print(aucs_tlr.shape)
    #mid_aucs = np.mean(aucs,axis=1)
    mid_aucs_tlr = np.mean(aucs_tlr,axis=1)
    mid_asdsfs = [np.mean(a) for a in asdsfs]
    print(mid_aucs_tlr.shape)


    # if dir=="primates_4_1000":  
    #     mid_aucs_tlr = mid_aucs_tlr[2:]
    #     mid_asdsfs = mid_asdsfs[2:]
    #     v = 1000
    #     for x, y in zip(mid_asdsfs,mid_aucs_tlr):
    #         v += 1000
    #         # ax.text(x*1.4, y, "%d" %v, ha="center")
    # if dir=="primates_1_200":
    #     mid_aucs_tlr = mid_aucs_tlr[:9]
    #     mid_asdsfs = mid_asdsfs[:9]
    #     v = 0
    #     for x, y in zip(mid_asdsfs,mid_aucs_tlr):
    #         v += 200
    #         # ax.text(x, y - 0.02, "%d" %v, ha="center")

    print(len(mid_asdsfs))
    print(len(mid_asdsfs))
    print(len(mid_aucs_tlr))

    iterations = np.loadtxt("primates_iterations.txt")
    inds = list(map(lambda x: int(x//200 -1),iterations))
    mid_asdsfs = [j for i,j in enumerate(mid_asdsfs) if i in inds]

    print(len(mid_asdsfs))
    plt.scatter(mid_asdsfs,mid_aucs_tlr,c=iterations,norm=matplotlib.colors.LogNorm())
    plt.colorbar()
    plt.xscale("log")
    plt.ylim([0.5,1.02])
    v = float(dir.split("_")[-1])
    
    #plt.plot(v * np.array(range(2,2+len(mid_aucs_tlr))),mid_aucs_tlr,"o",label=f"diagnfreq={int(v)}")
    #plt.semilogx(cs,mid_aucs_tlr,"o",label=f"diagnfreq={int(v)}")
    #plt.semilogy(v * np.array(range(2,2+len(mid_asdsfs))),mid_asdsfs,"o",label=dir)

    # for i in range(15):
    #     a = asdsfs[i,:]
    #     b = aucs[i,:]
    #     print(np.std(b))
    #     a = a[a>0]
        #plt.hist(b,bins=30)
    #     # plt.show()
    #plt.plot(5000*np.array(range(1,1+len(aucs))), aucs,"ok")
    #plt.semilogy(5000*np.array(range(1,1+len(asdsfs))), asdsfs,"ok")

# plt.legend()
plt.xlabel("ASDSF")
plt.ylabel("AUC")
plt.show()




# for i in range(N_aucs):
#     with open(os.path.join(gene_folder,"aucs", f'auc_{i}.npy'), 'rb') as f:
#         ngens = np.load(f)
#         aucvals = np.load(f)
#         aucs[i] = aucvals
# for i in range(N_asdsfs):
#     iters = np.loadtxt(os.path.join(gene_folder,"iters",f"iter_{i}"))
#     asdsf = [np.loadtxt(os.path.join(gene_folder,"asdsfs",f"asdsf_{i}"))[j] for j in sorted(range(len(iters)), key= lambda x: iters[x])]
#     asdsfs[i] = asdsf
# aucs = np.max(aucs,axis=0)
# asdsfs = np.max(asdsfs,axis=0)
# # q=1
# # aucs = np.quantile(aucs,q,axis=0)
# # asdsfs = np.quantile(asdsfs,q,axis=0)
# # print(aucs)
# # print(asdsfs)
# all_aucs[k] = aucs
# all_asdsfs[k] = asdsfs

# # plt.semilogx(ngens,aucs)
# # plt.semilogx(ngens,asdsfs)
# # plt.xlabel("Number of generations")
# # plt.show()

# mx = np.min(all_asdsfs)
# Mx = np.max(all_asdsfs)
# my = np.min(all_aucs)
# My = np.max(all_aucs)

# fair_excellent_threshold_asdsfs= 0.05
# bad_fair_threshold_asdsfs= 0.01 
# bad_fair_threshold_aucs = 0.6
# fair_excellent_threshold_aucs = 0.8

# mat = np.zeros((3,3))

# asdf_thresholds = [-float("inf"),bad_fair_threshold_asdsfs,fair_excellent_threshold_asdsfs,float("inf")]
# auc_thresholds = [-float("inf"),bad_fair_threshold_aucs,fair_excellent_threshold_aucs,float("inf")]
# for i in range(3):
#     l_asdfs = np.logical_and(all_asdsfs>=asdf_thresholds[i],all_asdsfs<asdf_thresholds[i+1] )
#     for j in range(3):
#         l_aucs = np.logical_and(all_aucs>=auc_thresholds[j],all_aucs<auc_thresholds[j+1] )
#         mat[-j-1][i]=int(sum(sum(np.logical_and(l_asdfs,l_aucs))))

# print(mat)


# # ngens=[10,20,30,40,50,75,100,150,200,300,500]
# # for j in range(m):    
# #     fig = plt.figure()
# #     ax = fig.add_subplot(111)
# #     asdsfs = all_asdsfs[:,j]
# #     aucs   = all_aucs[:,j]
# #     # for x, y, j in zip(asdsfs,aucs,range(num_genes)):
# #     #     v = ngens[j]
# #     #     ax.text(x * 1.2, y, "%d" %v, ha="center")

# #     plt.semilogx(asdsfs,aucs,"ok")
# #     plt.plot([bad_fair_threshold_asdsfs,bad_fair_threshold_asdsfs],[my,My],"--k")
# #     plt.plot([fair_excellent_threshold_asdsfs,fair_excellent_threshold_asdsfs],[my,My],"--k")
# #     plt.plot([mx,Mx],[bad_fair_threshold_aucs,bad_fair_threshold_aucs],"--k")
# #     plt.plot([mx,Mx],[fair_excellent_threshold_aucs,fair_excellent_threshold_aucs],"--k")
# #     plt.xlabel("ASDSF")
# #     plt.ylabel("AUC")
# #     plt.title(f"Ngen={ngens[j]}")
# #     plt.show(block=False)
# #     plt.pause(1)
# #     plt.close()

# # for i in range(num_genes):    
# #     fig = plt.figure()
# #     ax = fig.add_subplot(111)
# #     asdsfs = all_asdsfs[i,:]
# #     aucs   = all_aucs[i,:]
# #     for x, y, j in zip(asdsfs,aucs,range(num_genes)):
# #         v = ngens[j]
# #         ax.text(x * 1.2, y, "%d" %v, ha="center")

# #     plt.semilogx(asdsfs,aucs,"ok")
# #     plt.plot([bad_fair_threshold_asdsfs,bad_fair_threshold_asdsfs],[my,My],"--k")
# #     plt.plot([fair_excellent_threshold_asdsfs,fair_excellent_threshold_asdsfs],[my,My],"--k")
# #     plt.plot([mx,Mx],[bad_fair_threshold_aucs,bad_fair_threshold_aucs],"--k")
# #     plt.plot([mx,Mx],[fair_excellent_threshold_aucs,fair_excellent_threshold_aucs],"--k")
# #     plt.xlabel("ASDSF")
# #     plt.ylabel("AUC")
# #     plt.title(f"Gene {i+1}")
# #     if i in [2,14]:
# #         plt.show()
# #     else:
# #         plt.show(block=False)
# #         plt.pause(1)
# #         plt.close()
    
# fig = plt.figure()
# ax = fig.add_subplot(111)
# for i in range(num_genes):    
#     asdsfs = all_asdsfs[i,:]
#     aucs   = all_aucs[i,:]
#     plt.semilogx(asdsfs,aucs,"o")
#     plt.plot([bad_fair_threshold_asdsfs,bad_fair_threshold_asdsfs],[my,My],"--k")
#     plt.plot([fair_excellent_threshold_asdsfs,fair_excellent_threshold_asdsfs],[my,My],"--k")
#     plt.plot([mx,Mx],[bad_fair_threshold_aucs,bad_fair_threshold_aucs],"--k")
#     plt.plot([mx,Mx],[fair_excellent_threshold_aucs,fair_excellent_threshold_aucs],"--k")
# plt.xlabel("ASDSF")
# plt.ylabel("AUC")
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# for j in range(m):    
#     asdsfs = all_asdsfs[:,j]
#     aucs   = all_aucs[:,j]
#     plt.semilogx(asdsfs,aucs,"o")
#     plt.plot([bad_fair_threshold_asdsfs,bad_fair_threshold_asdsfs],[my,My],"--k")
#     plt.plot([fair_excellent_threshold_asdsfs,fair_excellent_threshold_asdsfs],[my,My],"--k")
#     plt.plot([mx,Mx],[bad_fair_threshold_aucs,bad_fair_threshold_aucs],"--k")
#     plt.plot([mx,Mx],[fair_excellent_threshold_aucs,fair_excellent_threshold_aucs],"--k")
# plt.xlabel("ASDSF")
# plt.ylabel("AUC")
# plt.show()

# for j in range(m):
#     plt.semilogy([ngens[j]]*len(all_asdsfs),all_asdsfs[:,j],"ok")
# plt.xlabel("Ngen")
# plt.title("ASDSF")
# plt.show()
# for j in range(m):
#     plt.plot([ngens[j]]*len(all_aucs),all_aucs[:,j],"ok")
# plt.xlabel("Ngen")
# plt.title("AUC")
# plt.show()