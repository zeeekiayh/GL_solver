#!/usr/bin/python3
# import math as m
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# from numpy.lib.function_base import copy
fig = plt.figure(figsize=(39,30))
gs = gridspec.GridSpec(3,4)

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
      # transpose the matrix and mirror across x and y to arrange it nicely to plot
      new_arr[int(size_x)-col-1][int(size_y)-row-1] = a # add the element to the 2D array
      
   return new_arr
      
data = np.loadtxt('data.txt')

# axes ranges/extent for imshow()
ext = [-size_x*step/2,size_x*step/2,0,size_y*step]

# convert all the OP components to 2D arrays
# and plot them in a grid of plots
ax1 = fig.add_subplot(gs[0,0])
Axx_real = unFlatten(data[:,2])
im = plt.imshow(Axx_real,vmin=-1,vmax=1,extent=ext)
ax1.set_title("Axx_real")

ax2 = fig.add_subplot(gs[0,1])
Axx_comp = unFlatten(data[:,3])
plt.imshow(Axx_comp,vmin=-1,vmax=1,extent=ext)
ax2.set_title("Axx_comp")

ax3 = fig.add_subplot(gs[0,2])
Axz_real = unFlatten(data[:,4])
plt.imshow(Axz_real,vmin=-1,vmax=1,extent=ext)
ax3.set_title("Axz_real")

ax4 = fig.add_subplot(gs[0,3])
Axz_comp = unFlatten(data[:,5])
plt.imshow(Axz_comp,vmin=-1,vmax=1,extent=ext)
ax4.set_title("Axz_comp")

ax5 = fig.add_subplot(gs[1,0])
Ayy_real = unFlatten(data[:,6])
plt.imshow(Ayy_real,vmin=-1,vmax=1,extent=ext)
ax5.set_title("Ayy_real")

ax6 = fig.add_subplot(gs[1,1])
Ayy_comp = unFlatten(data[:,7])
plt.imshow(Ayy_comp,vmin=-1,vmax=1,extent=ext)
ax6.set_title("Ayy_comp")

ax7 = fig.add_subplot(gs[1,2])
Azx_real = unFlatten(data[:,8])
plt.imshow(Azx_real,vmin=-1,vmax=1,extent=ext)
ax7.set_title("Azx_real")

ax8 = fig.add_subplot(gs[1,3])
Azx_comp = unFlatten(data[:,9])
plt.imshow(Azx_comp,vmin=-1,vmax=1,extent=ext)
ax8.set_title("Azx_comp")

ax9 = fig.add_subplot(gs[2,0])
Azz_real = unFlatten(data[:,10])
plt.imshow(Azz_real,vmin=-1,vmax=1,extent=ext)
ax9.set_title("Azz_real")

ax10 = fig.add_subplot(gs[2,1])
Azz_comp = unFlatten(data[:,11])
plt.imshow(Azz_comp,vmin=-1,vmax=1,extent=ext)
ax10.set_title("Azz_comp")

plt.xlabel('x')
plt.ylabel('z')

# colorbar placement and size:
#                       x_pos, y_pos width, height
cbar_ax = fig.add_axes([0.85,  0.1, 0.02,  0.25])
fig.colorbar(im, cax=cbar_ax)

plt.savefig('He3Defect.png')
# plt.show()