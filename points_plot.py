#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(16,8))
gs = gridspec.GridSpec(2,2)

conditions = []
for line in open('conditions.txt'):
   if not line.startswith('#'):
      try: conditions.append(float(line))
      except ValueError: continue # not a number, so disregard
x = np.linspace(0,conditions[0]*conditions[2],int(conditions[0]))

ax = fig.add_subplot(gs[:,0])

points_real1 = np.loadtxt('para_real.txt')
plt.plot(x,points_real1,label='$Re(\Delta_{\parallel})$')

points_comp1 = np.loadtxt('para_comp.txt')
plt.plot(x,points_comp1,label='$Im(\Delta_{\parallel})$')

points_real2 = np.loadtxt('perp_real.txt')
plt.plot(x,points_real2,label='$Re(\Delta_{\perp})$')

points_comp2 = np.loadtxt('perp_comp.txt')
plt.plot(x,points_comp2,label='$Im(\Delta_{\perp})$')

# mark the important points
# ax.annotate('$\sqrt{3}$',
#             xy=(3**(1/2),-0.01),
#             xytext=(3**(1/2),-0.15),
#             arrowprops=dict(facecolor='black',shrink=0.1))
plt.title('Order Parameter')
plt.xlabel(r'$z/\xi_0$')
plt.ylabel('$\Delta/\Delta_0$')
plt.legend()

ax = fig.add_subplot(gs[0,1])
integ_real = np.loadtxt('integ_r.txt')
integ_comp = np.loadtxt('integ_c.txt')
plt.plot(np.linspace(0,conditions[0]*conditions[2],len(integ_real)),integ_real,label='real')
plt.plot(np.linspace(0,conditions[0]*conditions[2],len(integ_comp)),integ_comp,label='imaginary')
plt.title('Integrand')
plt.legend()

ax = fig.add_subplot(gs[1,1])
plt.title('Complex plane')
plt.plot(points_real1,points_comp1,label='$\Delta_1$')
plt.plot(points_real2,points_comp2,label='$\Delta_2$')
plt.xlabel('real')
plt.ylabel('imaginary')
plt.legend()

plt.show()