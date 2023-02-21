# d = [
#     [0,5,3,8],
#     [5,0,6,1],
#     [3,6,0,2],
#     [8,1,2,0]
# ]

# static class DisjointSet {
# 		int[] s;
		
# 		public DisjointSet(int n) {
# 				Arrays.fill(s = new int[n], -1);
# 		}

# 		public int find(int i) {
# 				return s[i] < 0 ? i : (s[i] = find(s[i]));
# 		}

# 		public boolean union(int a, int b) {
# 				if ((a = find(a)) == (b = find(b))) return false;
# 				if (s[a] == s[b]) s[a]--;
# 				if (s[a] <= s[b]) s[b] = a; else s[a] = b;
# 				return true;
# 		}
# }
# def find(x):
#     if x != parents[x]:
#         parents[x] = find(parents[x])
#     return parents[x]    
# def union(x, y):
#     parents[find(x)] = find(y)

# import numpy as np
def um_projection(x):
    e = len(x)
    m = (1+(1+8*e)**0.5)//2 # solution of the quadratic (m choose 2) = e
    m = int(m)
    mat_to_vec = lambda i,j:    i * (2*m-i-1)//2 + j - i -1 if j > i else mat_to_vec(j,i)
    def vec_to_mat(a):
        i = int( (m-0.5) - ( (m-0.5)**2 - 2 * a )**0.5 ) 
        j = a - i * m + (i+1)*(i+2)//2
        return i,j 
    
    def find(x):
        if s[x] < 0:    return x
        s[x] = find(s[x])
        return s[x]

    def union(a,b,v): # returns the number of elements in the new union class
        a,b = find(a),find(b)
        if (a==b):  return len(d[a])
        if s[a] == s[b]:    s[a] -= 1
        if s[a] > s[b]: b,a=a,b 
        s[b] = a 
        for i in d[a]:
            for j in d[b]: 
                proj[mat_to_vec(i,j)] = v
        d[a] += d[b]   
        return len(d[a])

    data = sorted([ (di,vec_to_mat(i)) for i,di in enumerate(x) ])
    s, d = [-1] * m, [ [i] for i in range(m) ]
    proj = [None] * e
    for v,(i,j) in data:
        if union(i,j,v) == m: return proj #np.array(proj) + 1/2 * max(abs(np.array(proj) - np.array(x) )) 

import time
# tic=time.perf_counter()
# d = [5,3,8,6,1,2]
# print(um_projection(d))
# toc= time.perf_counter()
# print(toc-tic)
# tic=time.perf_counter()
# print(max(d))
# toc= time.perf_counter()
# print(toc-tic)

# d = [17,21,31,23,30,34,21,28,39,43]
d = [5,3,8,6,1,2]
tic=time.perf_counter()
um_projection(d)
toc= time.perf_counter()
print(toc-tic)
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
# dn = dendrogram(linkage(d))
# plt.show()
tic=time.perf_counter()
linkage(d)
toc= time.perf_counter()
print(toc-tic)