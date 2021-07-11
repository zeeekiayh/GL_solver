#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# fig = plt.figure(figsize=(16,8))
# gs = gridspec.GridSpec(2,2)

conditions = []
for line in open('conditions.txt'):
   if not line.startswith('#'):
      try: # CONTINUE HERE
         # NEED TO CHANGE THE WAY i READ IN THE VALUES
         num = float(line)
         print(num)
         conditions.append(num)
      except ValueError: continue # not a number, so disregard
# x = np.linspace(0,conditions[0]*conditions[2],int(conditions[0]))

def unFlatten(arr):
   # conditions[2]: size of the mesh (on a side)
   copy_arr = []
   counts = 0
   sub_list = [] # hold the sub lists to add to arr
   print(conditions[2])
   for a in arr:
      print('\t'+str(a)+'\tcounts = '+str(counts))
      sub_list.append(a)
      if (counts+1)%conditions[2] == 0:
         copy_arr.append([sub_list])
         sub_list = []
      counts += 1
   
   return copy_arr
      
# ax = fig.add_subplot(gs[:,0])
data = np.loadtxt('data.txt')

Axx_real = data[:,2]
Axx_comp = data[:,3]

Axz_real = data[:,4]
Axz_comp = data[:,5]

Ayy_real = data[:,6]
Ayy_comp = data[:,7]

Azx_real = data[:,8]
Azx_comp = data[:,9]

Azz_real = data[:,10]
Azz_comp = data[:,11]

print(Axx_real)
print(unFlatten(Axx_real))

plt.show()