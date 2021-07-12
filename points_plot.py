#!/usr/bin/python3
import math as m
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.lib.function_base import copy
# fig = plt.figure(figsize=(16,8))
# gs = gridspec.GridSpec(2,2)

# read in the conditions from the file
conditions = []
for line in open('conditions.txt'):
   if not line.startswith('#'): # lines starting with '#' are just titles
      try:
         num = list(map(float, line.split()))  # put all the values into a list
         if len(num)>0: conditions.append(num) # only add if the list is not empty
      except ValueError: continue # not a number, so disregard

# define variables from the conditions list
size_x = conditions[0][0]
size_y = conditions[0][1]
step   = conditions[2][0]

# a function to convert a 1D array to a 2D
def unFlatten(arr):
   new_arr = np.zeros((int(size_x),int(size_y))) # the list to return
   row = 0  # row number
   col = -1 # column number

   for a in arr: # look at all elements in the array
      if (col+1)%size_x == 0 and not col == -1: # if we're at the last element...
         row += 1 # we go to the next row
         col = -1 # and reset the column number
      col += 1
      new_arr[row][col] = a # add the element to the 2D array
      
   return new_arr
      
# ax = fig.add_subplot(gs[:,0])
data = np.loadtxt('data.txt')

# convert all the OP components to 2D arrays
Axx_real = unFlatten(data[:,2])
Axx_comp = unFlatten(data[:,3])

Axz_real = unFlatten(data[:,4])
Axz_comp = unFlatten(data[:,5])

Ayy_real = unFlatten(data[:,6])
Ayy_comp = unFlatten(data[:,7])

Azx_real = unFlatten(data[:,8])
Azx_comp = unFlatten(data[:,9])

Azz_real = unFlatten(data[:,10])
Azz_comp = unFlatten(data[:,11])

plt.imshow(Axx_real)
plt.show()