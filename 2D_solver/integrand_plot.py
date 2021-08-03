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

plt.plot(np.linspace(0,step*size_z,len(data)),data)
plt.show()