import matplotlib.pyplot as plt
import numpy as np
import time as tm
import os

# CODE FOR RUNNING MANY SYSTEMS FOR DIFFERENT r
r_start = 8.5
r_step = 0.5
r_end = 40

start_time_seconds = tm.time()
start_time_pretty = tm.ctime()
print(f'Start time: {start_time_pretty}')

# for radius in np.arange(r_start, r_end, r_step):
radius = r_start
while radius <= r_end:
    print(f'\n====================================\nRunning "gl_fdm" for {radius = }')
    os.system(f'./gl_fdm conditions5c_bubble.txt c {radius} y')
    radius += r_step
    if radius > 17.9:
        r_step = 1
print(f'\n====================================\n')

end_time_seconds = tm.time()
end_time_pretty = tm.ctime()
print(f'Total runtime: {end_time_seconds - start_time_seconds}')
print(f'Start time: {start_time_pretty}')
print(f'End time: {end_time_pretty}')


# CODE FOR PLOTTING THE RESULTS
file_name1 = 'E_vs_r.txt'
data1 = np.loadtxt(file_name1)

file_name2 = 'E_vs_r_size20.txt'
data2 = np.loadtxt(file_name2)

# check to see if there are duplicates?
# sort the list at all? so normal plotting looks good

power=1
# plt.plot(data1[:,1]**power,data1[:,0],label='data1')
# plt.plot(data2[:,1]**power,data2[:,0],label='data2')
plt.scatter(data2[:,1]**power,data2[:,0],label='for h=0.5')
plt.scatter(data1[:,1]**power,data1[:,0],label='for h=0.13')
# plt.xlim(0,8)
# plt.ylim(0,250)
plt.legend()
plt.xlabel(rf'wall radius ($(r/\xi_0)$)')#^{power}
plt.ylabel('Energy')
plt.title('energy of OP configuration relative to refFE')
plt.show()