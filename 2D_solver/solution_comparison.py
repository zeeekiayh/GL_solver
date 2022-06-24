import sys # for arguments in terminal call
import numpy as np
import matplotlib.pyplot as plt

# to return "Nop, N_FE_cols, Nu, Nv, h, X, Z, OP_data_array, FE_data_array, labels"
def read_file(file_name):
    OP_data_array = [] # hold all OP components; to be converted to np.array() later
    FE_data_array = [] # holds all FE data
    X, Z = [], []      # the positions arrays
    labels = []        # to store the labels of all the columns

    # counters for grid size
    u_0_count = 0
    v_0_count = 0

    # initialize Nop and N_FE_cols
    Nop = 0
    N_FE_cols = 0

    line_count = 0
    for line in open(file_name,'r'):
        # the first line should give us all the labels of the columns
        if line_count == 0:
            labels = line.split()
            # loop through all the labels and determine Nop
            for l in labels:
                if l[0] == "#":
                    Nop += 1
            Nop //= 2 # to account for complex columns

            # number of columns excluding position and OP component columns
            N_FE_cols = len(labels) - 2 - 2*Nop
        else:
            # get all the values from the line
            values = list(map(float,line.split()))

            # the first 2 values are for X and Z
            X.append(values[0])
            Z.append(values[1])

            # start reading in the values into OP_data_array
            OP_data_array.append( values[2:-N_FE_cols] )

            # read in the values into FE_data_array
            FE_data_array.append( values[-N_FE_cols:] )

            # make counts to determine grid dimensions
            if (abs(X[-1]) < 1e-8):
                u_0_count += 1
            if (abs(Z[-1]) < 1e-8):
                v_0_count += 1
        line_count += 1

    Nu = int(v_0_count) # YES! the u and v are switched! they need to be!
    Nv = int(u_0_count)

    # convert the data array and reshape all components
    OP_data_array = np.array(OP_data_array)
    OP_data_array = np.array([ np.reshape(OP_data_array[:,2*i], (Nv,Nu)).transpose() for i in range(Nop) ])

    # step size; we'll assume it will be the same in all directions
    h = X[1]-X[0]
    # in case position values are in a different order...
    if abs(h) < 1e-8: h = Z[1]-Z[0]

    # reshape the X and Z arrays
    X = np.reshape(np.array(X), (Nv,Nu))
    Z = np.reshape(np.array(Z), (Nv,Nu))

    # convert the data array
    FE_data_array = np.array(FE_data_array)
    FE_data_array = np.array([ np.reshape(FE_data_array[:,i], (Nv,Nu)).transpose() for i in range(N_FE_cols) ])
        
    return Nop, N_FE_cols, Nu, Nv, h, X, Z, OP_data_array, FE_data_array, labels

# plot all the OP components in 2D and slices
def plot_OP_comps_and_slices_compare(file_name1, file_name2):
    Nop1, N_FE_cols1, Nu1, Nv1, h1, X1, Z1, OP_data_array1, FE_data_array1, labels1 = read_file(file_name1)
    Nop2, N_FE_cols2, Nu2, Nv2, h2, X2, Z2, OP_data_array2, FE_data_array2, labels2 = read_file(file_name2)

    if input(f'Is {file_name1} for a cylindrical system? (y/n): ') == 'y': custom_labels1 = [r'$A_{rr}$', r'$A_{\phi \phi}$', r'$A_{zz}$', r'$A_{zr}$', r'$A_{rz}$'] # cylindrical
    else: custom_labels1 = [r'$A_{xx}$', r'$A_{yy}$', r'$A_{zz}$', r'$A_{zx}$', r'$A_{xz}$']

    if input(f'Is {file_name2} for a cylindrical system? (y/n): ') == 'y': custom_labels2 = [r'$A_{rr}$', r'$A_{\phi \phi}$', r'$A_{zz}$', r'$A_{zr}$', r'$A_{rz}$'] # cylindrical
    else: custom_labels2 = [r'$A_{xx}$', r'$A_{yy}$', r'$A_{zz}$', r'$A_{zx}$', r'$A_{xz}$']

    # just assume that these are all the same
    Nop, N_FE_cols, Nu, Nv, h, X, Z = Nop1, N_FE_cols1, Nu1, Nv1, h1, X1, Z1

    # But check to see if the h's are the same! That would be a difficult one to see in the plots!
    if (abs(h1-h2) > 1e-8): print(f'WARNING: Your step sizes are noticeably different: {h1 - h2 = }; this could cause larger differences than expected!')

    try:
        # calculate the difference between the arrays (assuming that everything else is the same)
        OP_data_array, FE_data_array = OP_data_array1-OP_data_array2, FE_data_array1-FE_data_array2
    except ValueError:
        print("ERROR: you probably gave files with differently shaped domains. Check your conditions files.")
        print(f'\tShape from {file_name1}: {np.shape(OP_data_array1)}')
        print(f'\tShape from {file_name2}: {np.shape(OP_data_array2)}')
        exit()

    # the domain extents for the imshow calls
    ext = [0*h, Nu*h, 0*h, Nv*h]

    # initialize all axes
    axs = []
    for _ in range(Nop): axs.append(None)
    # other axes we may use ... do we need to add more?
    FE_ax, FE_prof_ax, empty_ax, grad_FE_ax = None, None, None, None
    fig, axes = None, None

    # shape the plot based on OP size
    if Nop == 3:
        fig, axes = plt.subplots(3,2)
        # then unpack the axes tuple
        ((grad_FE_ax, axs[0]),
         (FE_ax,      axs[1]),
         (FE_prof_ax, axs[2])) = axes

        if Nu > 1 and Nv > 1:
            axs[2].set_xlabel(rf'$x/\xi_0$ (bottom/surface)')
            axs[0].axes.xaxis.set_ticks([])
            axs[1].axes.xaxis.set_ticks([])
            axs[0].axes.yaxis.set_ticks([])
            axs[1].axes.yaxis.set_ticks([])
        elif Nu == 1:
            axs[0].set_xlabel(rf'$z/\xi_0$')
            axs[1].set_xlabel(rf'$z/\xi_0$')
            axs[2].set_xlabel(rf'$z/\xi_0$')
        elif Nv == 1:
            axs[0].set_xlabel(rf'$x/\xi_0$')
            axs[1].set_xlabel(rf'$x/\xi_0$')
            axs[2].set_xlabel(rf'$x/\xi_0$')
    elif Nop == 5:
        fig, axes = plt.subplots(3,3)
        # then unpack the axes tuple
        ((grad_FE_ax, axs[0], axs[3]),
         (FE_ax,      axs[1], axs[4]),
         (FE_prof_ax, axs[2], empty_ax)) = axes

        axs[0].axes.xaxis.set_ticks([])
        axs[1].axes.xaxis.set_ticks([])
        axs[3].axes.xaxis.set_ticks([])

        if Nu > 1 and Nv > 1:
            axs[0].axes.yaxis.set_ticks([])
            axs[1].axes.yaxis.set_ticks([])
            axs[3].axes.yaxis.set_ticks([])
            axs[4].axes.yaxis.set_ticks([])

            axs[2].set_xlabel(x_axis_labels)
            axs[4].set_xlabel(x_axis_labels)
        elif Nu == 1:
            axs[2].set_xlabel(rf'$z/\xi_0$')
            axs[4].set_xlabel(rf'$z/\xi_0$')
        elif Nv == 1:
            axs[2].set_xlabel(rf'$x/\xi_0$')
            axs[4].set_xlabel(rf'$x/\xi_0$')
    else: print(f"Implement 'plot_OP_comps_and_slices' for {Nop = }.")

    fig.suptitle(f'OP-{Nop}: Diff btwn {file_name1} & {file_name2}')

    # plot the 2D solution
    for i in range(Nop):
        axs[i].set_title(f'{custom_labels1[i]} - {custom_labels2[i]}')
        if Nu > 1 and Nv > 1:
            im = axs[i].imshow(OP_data_array[i][:,::-1].transpose(), extent=ext, cmap='bwr')
            plt.colorbar(im,ax=axs[i])
        elif Nu == 1:
            axs[i].plot(np.linspace(ext[2],ext[3],Nv), OP_data_array[i][0])
            axs[i].set_ylim( [ min([0, OP_data_array.min()])-0.1, 0.1+max([0, OP_data_array.max()]) ] )
        elif Nv == 1:
            axs[i].plot(np.linspace(ext[0],ext[1],Nu), OP_data_array[i])
            axs[i].set_ylim( [ min([0, OP_data_array.min()])-0.1, 0.1+max([0, OP_data_array.max()]) ] )
    
    if empty_ax != None:
        empty_ax.axis('off')

    # plot the FE profile
    FE_prof_ax.set_ylabel('FE')
    x_array = np.linspace(ext[2],ext[3],Nv)
    if Nu > 1 and Nv > 1:
        FE_prof_ax.set_title(rf'Total FE profile; along $x/\xi_0 = {round(x_array[len(FE_data_array[0])//2],2)}$')
        FE_prof_ax.plot(x_array, FE_data_array[0][len(FE_data_array[0])//2,:])
        FE_prof_ax.set_xlabel(r'$z/\xi_0$')
    elif Nu == 1:
        FE_prof_ax.set_title(r'Total FE profile; along $x/\xi_0 = 0$')
        FE_prof_ax.plot(x_array, FE_data_array[0][0])
        FE_prof_ax.set_xlabel(rf'$z/\xi_0$')
    elif Nv == 1:
        x_array = np.linspace(ext[0],ext[1],Nv)
        FE_prof_ax.set_title(r'Total FE profile; along $z/\xi_0 = 0$')
        FE_prof_ax.plot(x_array, FE_data_array[0])
        FE_prof_ax.set_xlabel(r'$x/\xi_0$')

    # plot the 2D heatmap of the FE
    FE_ax.set_title('Total FE')
    FE_ax.axes.xaxis.set_ticks([])
    if Nu > 1 and Nv > 1:
        FE_ax.set_ylabel(rf'$z/\xi_0$ (left)')
        im = FE_ax.imshow(FE_data_array[0][:,::-1].transpose(), extent=ext, cmap='gist_heat')
        fig.colorbar(im,ax=FE_ax)
    elif Nu == 1:
        FE_ax.set_ylabel('amplitude')
        FE_ax.plot(np.linspace(ext[2],ext[3],Nv), FE_data_array[0][0])
    elif Nv == 1:
        FE_ax.set_ylabel('amplitude')
        FE_ax.plot(np.linspace(ext[0],ext[1],Nu), FE_data_array[0])

    # plot the defect energy
    grad_FE_ax.set_title('Grad free energy')
    grad_FE_ax.axes.xaxis.set_ticks([])
    if Nu > 1 and Nv > 1:
        grad_FE_ax.set_ylabel(rf'$z/\xi_0$ (left)')
        im = grad_FE_ax.imshow(FE_data_array[2][:,::-1].transpose(), extent=ext, cmap='gist_heat')
        fig.colorbar(im,ax=grad_FE_ax)
    elif Nu == 1:
        grad_FE_ax.set_ylabel('amplitude')
        grad_FE_ax.plot(np.linspace(ext[2],ext[3],Nv), FE_data_array[2][0])
    elif Nv == 1:
        grad_FE_ax.set_ylabel('amplitude')
        grad_FE_ax.plot(np.linspace(ext[0],ext[1],Nu), FE_data_array[2])

    # display all plots plotted above
    plt.show()

# the main code to run
def main(argv):
    # argv should be like: [ file_name1, file_name2 ]
    plot_OP_comps_and_slices_compare(argv[0], argv[1]) # pass file names to function

# how the python interpreter will know to run our function called 'main()'
if __name__ == "__main__":
    main(sys.argv[1:]) # call the main function with the arguments given in the terminal call
