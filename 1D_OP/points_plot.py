#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

conditions = np.loadtxt('conditions.txt')
x = np.linspace(0,conditions[0]*conditions[2],int(conditions[0]))

try:
   points = np.loadtxt('points.txt')
   plt.plot(x,points)
except:
   print('points.txt file not found.')
# plt.legend()
plt.show()