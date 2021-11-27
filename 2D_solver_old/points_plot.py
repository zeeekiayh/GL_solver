#!/usr/bin/python3
# import math as m
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# from numpy.lib.function_base import copy
fig = plt.figure(figsize=(30,20))
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
   new_arr = np.zeros((int(size_y),int(size_x))) # the list to return
   row_number = 0  # row number
   col_number = -1 # column number

   for a in arr: # look at all elements in the array
      if (col_number+1)%size_x == 0 and not col_number == -1: # if we're at the last element in a row
         row_number += 1 # we go to the next row
         col_number = -1 # and reset the column number
      col_number += 1
      if row_number < size_y and col_number < size_x:
         new_arr[row_number][col_number] = a # add the element to the 2D array
      
   return new_arr
      
data = np.loadtxt('data.txt')
# data = np.loadtxt('test_solution.txt')

# axes ranges/extent for imshow()
ext = [min(data[:,0]),max(data[:,0]),max(data[:,1]),min(data[:,1])]
x_range = np.linspace(0,size_y*step,int(size_y))
m.rc('font',size=20)

# convert all the OP components to 2D arrays
# and plot them in a grid of plots
ax1 = fig.add_subplot(gs[0,3])
Axx_real = unFlatten(data[:,2])
im1=plt.imshow(Axx_real,vmin=0,extent=ext)
plt.xlabel(r'$x/\xi_{||}$ (top)')
plt.ylabel(r'$z/\xi_\perp$ (right)')
ax1.set_title("Axx",x=0.9,y=1)
cbar_ax = fig.add_axes([0.9,0.65,0.02,0.22])
fig.colorbar(im1,cax=cbar_ax)

ax2 = fig.add_subplot(gs[1,3])
Ayy_real = unFlatten(data[:,3])
im2=plt.imshow(Ayy_real,vmin=0,extent=ext)
plt.xlabel(r'$x/\xi_{||}$ (top)')
plt.ylabel(r'$z/\xi_\perp$ (right)')
ax2.set_title("Ayy",x=0.9,y=1)
cbar_ax = fig.add_axes([0.9,0.39,0.02,0.22])
fig.colorbar(im2,cax=cbar_ax)

ax3 = fig.add_subplot(gs[2,3])
Azz_real = unFlatten(data[:,4])
im3=plt.imshow(Azz_real,vmin=0,extent=ext)
plt.xlabel(r'$x/\xi_{||}$ (top)')
plt.ylabel(r'$z/\xi_\perp$ (right)')
ax3.set_title("Azz",x=0.9,y=1)

# colorbar placement and size:
#                       x_pos, y_pos width, height
cbar_ax = fig.add_axes([0.9,  0.12,  0.02,  0.22])
fig.colorbar(im3,cax=cbar_ax)

# plot a slice of the op on the mesh
ax4 = fig.add_subplot(gs[:,:3])
plt.plot(x_range,Axx_real[:int(size_y),int(size_y/2)],'r',label='Axx')
plt.plot(x_range,Ayy_real[:int(size_y),int(size_y/2)],'-.b',label='Ayy')
plt.plot(x_range,Azz_real[:int(size_y),int(size_y/2)],'y',label='Azz')
plt.xlabel(r'$z/\xi_\perp$')
plt.ylabel(r'$\Delta/\Delta_0$')
plt.legend()

plt.savefig('He3Defect.png')
plt.show()