import matplotlib.pyplot as plt
import numpy as np
import time as tm
import os

# # CODE FOR RUNNING MANY SYSTEMS FOR DIFFERENT r
# r_start = 0.8
# r_step = 0.1
# r_end = 8.01

# start_time_seconds = tm.time()
# start_time_pretty = tm.ctime()
# print(f'Start time: {start_time_pretty}')

# for radius in np.arange(r_start, r_end, r_step):
#     print(f'\n====================================\nRunning "gl_fdm" for {radius = }')
#     os.system(f'./gl_fdm conditions5c_bubble.txt c {radius} y')
# print(f'\n====================================\n')
# end_time_seconds = tm.time()
# end_time_pretty = tm.ctime()
# print(f'Total runtime: {end_time_seconds - start_time_seconds}')
# print(f'Start time: {start_time_pretty}')
# print(f'End time: {end_time_pretty}')



# CODE FOR PLOTTING THE RESULTS
file_name = 'E_vs_r.txt'
data = np.loadtxt(file_name)

# check to see if there are duplicates?
# sort the list at all? so normal plotting looks good

power=1
plt.plot(data[:,1]**power,data[:,0])
# plt.scatter(data[:,1]**power,data[:,0])
plt.xlabel(rf'wall radius ($(r/\xi_0)^{power}$)')
plt.ylabel('Energy')
plt.title('energy of OP configuration relative to refFE')
plt.show()