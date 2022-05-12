#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(10,10))

# read in the conditions from the file
conditions = []
for line in open('conditions.txt'):
   if not line.startswith('#'): # lines starting with '#' are just titles
      try:
         num = list(map(float, line.split()))  # put all the values into a list
         if len(num)>0: conditions.append(num) # only add if the list is not empty
      except ValueError: continue # not a number, so disregard

# define variables from the conditions list
size_z = conditions[0][1]
step   = conditions[2][0]

data = np.loadtxt('integrand.txt')

plt.plot(data[:,0],data[:,1])
plt.title('total integrand')
plt.show()

integrand_bulk = np.loadtxt("integ_bulk.txt")
integrand_grad = np.loadtxt('integ_grad.txt')
integ_bulk_real = integrand_bulk[:,1]
integ_bulk_comp = integrand_bulk[:,2]
integ_grad_real = integrand_grad[:,1]
integ_grad_comp = integrand_grad[:,2]

plt.clf()
plt.plot(integ_bulk_real,label='bulk real')
plt.plot(integ_bulk_comp,label='bulk comp')
plt.plot(integ_grad_real,label='grad real')
plt.plot(integ_grad_comp,label='grad comp')
plt.legend()
plt.show()