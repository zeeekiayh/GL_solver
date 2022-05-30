import sys # for arguments in terminal call
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# returns "Nop, Nu, Nv, h, np_data_array, labels" from the solution file
def readSolutionFile(file_name):
    np_data_array = [] # to be converted to np.array() later
    labels = [] # to store the labels of all the columns

    u_0_count = 0
    v_0_count = 0

    line_count = 0
    for line in open(file_name,'r'):
        # the first line should give us all the labels of the columns
        if line_count == 0:
            labels = line.split()
        else:
            # start reading in the values into np_data_array
            np_data_array.append(list(map(float,line.split())))

            if (abs(np_data_array[-1][0]) < 1e-8):
                u_0_count += 1
            if (abs(np_data_array[-1][1]) < 1e-8):
                v_0_count += 1
        line_count += 1

    # convert the data array
    np_data_array = np.array(np_data_array)
    # this is now in the exact for that we had in the other python plotting file!

    # step size; we'll assume it will be the same in all directions
    h = np_data_array[1,0] - np_data_array[0,0]
    if abs(h) < 1e-8: h = np_data_array[1,1] - np_data_array[0,1] # in case position values are in a different order

    Nop = (len(labels)-2)//2 # take off the 2 cols for location; /2 b/c of complex cols

    Nu = int(v_0_count) # YES! the u and v are switched! they need to be!
    Nv = int(u_0_count)
    # print(f'{Nu = }; {Nv = }')
        
    return Nop, Nu, Nv, h, np_data_array, labels

# returns "Ncols, Nu, Nv, h, FE_data_array, labels" from the FE file
def read_FE_File(file_name):
    FE_data_array = [] # to be converted to np.array() later
    labels = [] # to store the labels of all the columns
    u_0_count = 0
    v_0_count = 0
    line_count = 0
    for line in open(file_name,'r'):
        # the first line should give us all the labels of the columns
        if line_count == 0:
            labels = line.split()
        else:
            # start reading in the values into FE_data_array
            FE_data_array.append(list(map(float,line.split())))
            if (abs(FE_data_array[-1][0]) < 1e-8):
                u_0_count += 1
            if (abs(FE_data_array[-1][1]) < 1e-8):
                v_0_count += 1
        line_count += 1

    # convert the data array
    FE_data_array = np.array(FE_data_array)
    # this is now in the exact for that we had in the other python plotting file!

    # step size; we'll assume it will be the same in all directions
    h = FE_data_array[1,0] - FE_data_array[0,0]
    if abs(h) < 1e-8: h = FE_data_array[1,1] - FE_data_array[0,1] # in case position values are in a different order

    # the number of columns in the FE file
    Ncols = int(len(labels)-2) # take off the 2 cols for location

    Nu = int(v_0_count) # YES! the u and v are switched! they need to be!
    Nv = int(u_0_count)
        
    return Ncols, Nu, Nv, h, FE_data_array, labels

# to return "Nop, N_FE_cols, Nu, Nv, h, X, Z, OP_data_array, FE_data_array, labels"
def read_file(file_name):
    # print("read_file")
    OP_data_array = [] # hold all OP components; to be converted to np.array() later
    FE_data_array = [] # holds all FE data
    X, Z = [], [] # the positions arrays
    labels = [] # to store the labels of all the columns

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

            # number of columns excluding position
            #   and OP component columns
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

    # print(labels[2:-N_FE_cols])
    # print(labels[-N_FE_cols:])
    # print(f'{Nu = }; {Nv = }; {Nop = }; {N_FE_cols = }')

    # convert the data array and reshape all components
    OP_data_array = np.array(OP_data_array)
    OP_data_array = np.array([ np.reshape(OP_data_array[:,2*i], (Nv,Nu)).transpose() for i in range(Nop) ])
    # OP_data_array = np.reshape( np.array(OP_data_array), (Nop*2,Nv,Nu) )

    # step size; we'll assume it will be the same in all directions
    h = X[1]-X[0]
    # in case position values are in a different order...
    if abs(h) < 1e-8: h = Z[1]-Z[0]

    # reshape the X and Z arrays
    X = np.reshape(np.array(X), (Nv,Nu))#.transpose()
    Z = np.reshape(np.array(Z), (Nv,Nu))#.transpose()

    # convert the data array
    FE_data_array = np.array(FE_data_array)
    FE_data_array = np.array([ np.reshape(FE_data_array[:,i], (Nv,Nu)).transpose() for i in range(N_FE_cols) ])
    # FE_data_array = np.reshape( np.array(FE_data_array), (N_FE_cols,Nv,Nu) )
        
    return Nop, N_FE_cols, Nu, Nv, h, X, Z, OP_data_array, FE_data_array, labels

# plot all the OP components in 2D and slices
def plot_OP_comps_and_slices(file_name):
    # print("plot_OP_comps_and_slices")
    Nop, N_FE_cols, Nu, Nv, h, X, Z, OP_data_array, FE_data_array, labels = read_file(file_name)
    custom_labels = [r'$A_{xx}$', r'$A_{yy}$', r'$A_{zz}$', r'$A_{zx}$', r'$A_{xz}$']

    # the domain extents for the imshow calls
    ext = [0*h, Nv*h, 0*h, Nu*h]

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

        axs[0].axes.xaxis.set_ticks([])
        axs[1].axes.xaxis.set_ticks([])
        axs[0].axes.yaxis.set_ticks([])
        axs[1].axes.yaxis.set_ticks([])
        axs[2].set_xlabel(rf'$z/\xi_0$ (right)')
    elif Nop == 5:
        fig, axes = plt.subplots(3,3)
        # then unpack the axes tuple
        ((grad_FE_ax, axs[0], axs[3]),
         (FE_ax,      axs[1], axs[4]),
         (FE_prof_ax, axs[2], empty_ax)) = axes

        axs[0].axes.xaxis.set_ticks([])
        axs[1].axes.xaxis.set_ticks([])
        axs[3].axes.xaxis.set_ticks([])
        axs[0].axes.yaxis.set_ticks([])
        axs[1].axes.yaxis.set_ticks([])
        axs[3].axes.yaxis.set_ticks([])
        axs[2].set_xlabel(rf'$z/\xi_0$ (right)')
        axs[4].set_xlabel(rf'$z/\xi_0$ (right)')
    else: print(f"Implement 'plot_OP_comps_and_slices' for {Nop = }.")

    fig.suptitle(f'OP-{Nop}')

    # plot the 2D solution
    for i in range(Nop):
        axs[i].set_title(f'{custom_labels[i]}')
        # axs[i].axes.xaxis.set_ticks([]) # remove tick marks for a cleaner plot;
        # axs[i].axes.yaxis.set_ticks([]) #   we only leave them on one for reference
        im = axs[i].imshow(OP_data_array[i], extent=ext, cmap='bwr') # .transpose()
        plt.colorbar(im,ax=axs[i])
    
    if empty_ax != None:
        empty_ax.axis('off') # don't show the blank plots
    #     # use the blank plot to show axes labels
    #     empty_ax.imshow(np.zeros(np.shape(OP_data_array[0])), extent=ext, vmin=-0.1, vmax=0.1, cmap='bwr')
    #     empty_ax.set_xlabel(rf'$z/\xi_0$ (right)')
    #     empty_ax.set_ylabel(rf'$x/\xi_0$ (bottom/surface)')
    #     empty_ax.set_title('axes labels')

    # plot the FE profile
    FE_prof_ax.set_xlabel(rf'$x/\xi_0$')
    FE_prof_ax.set_ylabel('FE')
    FE_prof_ax.set_title('Total FE profile')
    FE_prof_ax.plot(np.linspace(ext[2],ext[3],Nu), list(reversed(FE_data_array[0].transpose()[len(FE_data_array[0])//2,:]))) # .transpose()

    # plot the 2D heatmap of the FE
    FE_ax.set_title('Total FE')
    FE_ax.axes.xaxis.set_ticks([])
    # FE_ax.axes.yaxis.set_ticks([])
    FE_ax.set_ylabel(rf'$x/\xi_0$ (bottom/surface)')
    im = FE_ax.imshow(FE_data_array[0], extent=[0,ext[1],0,ext[3]], cmap='gist_heat') # .transpose()
    fig.colorbar(im,ax=FE_ax)

    # plot the defect energy
    grad_FE_ax.set_title('Grad free energy')
    grad_FE_ax.axes.xaxis.set_ticks([])
    # grad_FE_ax.axes.yaxis.set_ticks([])
    grad_FE_ax.set_ylabel(rf'$x/\xi_0$ (bottom/surface)')
    im = grad_FE_ax.imshow(FE_data_array[2], extent=[0,ext[1],0,ext[3]], cmap='gist_heat') # .transpose()
    fig.colorbar(im,ax=grad_FE_ax)

    # display all plots plotted above
    plt.show()

    # # FOR SOME REASON THIS 3D PLOT IS NOT WORKING ANY MORE...
    # # plot 3D, all components
    # plt.title(f'OP-{Nop}')
    # ax = plt.subplot(111, projection='3d')
    # ax.set_xlabel(r'$x/\xi_0$')
    # ax.set_ylabel(r'$z/\xi_0$')
    # ax.set_zlabel(r'$|A_{\alpha i}|$')
    # for i in range(Nop):
    #     # ax.scatter(X,Z,OP_data_array[i],label=f'OP comp {labels[i]}')
    #     # If the 3D plot is too crowded, use thinned out arrays!
    #     ax.scatter(X[::2,::2],Z[::2,::2],OP_data_array[i][::2,::2],label=f'OP comp #{labels[i]}')
    # plt.legend()
    # plt.show()

    # plt.clf()
    # # plot slices for 1D view
    # plt.clf()
    # plt.title("Slices top to bottom")
    # for i in range(Nop):
    #     print(f'{np.shape(organized_array) = }')
    #     slc = organized_array[i][:,len(organized_array[i][0])//2]
    #     plt.plot( np.linspace(0, h*size, len(slc)), slc, label=f"comp {i}" )
    # plt.legend()
    # plt.show()

# the main code to run
def main(argv):
    # argv should be like:
    #   [ Nop, [debug=false] ]
    debug=False
    if (len(argv)==2):
        print("Using debugging in python code...")
        debug=True
    elif (len(argv)!=1):
        print("Incorrect number of arguments in python call!")
        exit()

    plot_OP_comps_and_slices(argv[0])

    # # read values from conditions file
    # Nop, Nu, Nv, h, OP_array, labels = readSolutionFile(argv[0])

    # # sort the OP_array for ease of plotting
    # A = np.array([ np.reshape(OP_array[:,2*i+2], (Nv,Nu)).transpose() for i in range(Nop) ])
    # X, Z = np.reshape(OP_array[:,0], (Nv,Nu)), np.reshape(OP_array[:,1], (Nv,Nu))
    # X = X.transpose() # make sure they are good for display!
    # Z = Z.transpose()

    # # the domain extents for the imshow calls
    # ext = [0*h, Nv*h, 0*h, Nu*h]

    ## I don't have this working with the 'read_file()' function yet...
    # # debugging: visualize the initial guess
    # if debug:
    #     # Initial guess -- for debugging
    #     plt.title("initial guess")
    #     # get the guess from the saved file from the C++ code in gl_fdm.cpp
    #     Nop, Nu, Nv, h, op, labels = readSolutionFile(f'initGuess{Nop}.txt')
    #     initOP = np.array([ np.reshape(op[:,2*i+2], (Nv,Nu)) for i in range(Nop) ])

    #     if Nu > 1:
    #         plot_OP_comps_and_slices( Nop, initOP, X, Z, ext, list(map(lambda l: l+"_guess", labels[2::2])), None )
    #     elif Nu == 1: # basically for the 1D case
    #         # plot only the slices, since the 2D view is not useful here
    #         plt.title("Slices top to bottom")
    #         for i in range(Nop):
    #             slc = initOP[i][0]
    #             plt.plot( np.linspace(0,h*Nv,len(slc)), slc, label=f"initOP[{i}]" )
    #         plt.legend()
    #         plt.show()

    # if Nu > 1:
        # plot_OP_comps_and_slices(Nop, A, X, Z, ext, labels[2::2], f'totalFE{Nop}.txt' )

        # if debug:
        #     # Plot the bulk free energy
        #     N_FE_cols, Nu, Nv, h, FE_bulk, labels = read_FE_File(f'bulkRHS_FE{Nop}.txt')

        #     plt.title(f'Bulk Free Energy for OP-{Nop}')
        #     im = plt.imshow( np.reshape(FE_bulk[:,2], (Nv,Nu)), extent=ext, cmap='gist_heat' )
        #     plt.colorbar(im)
        #     plt.show()
            
        #     # Plot the grad free energy
        #     N_FE_cols, Nu, Nv, h, FE_grad, labels = read_FE_File(f'FEgrad{Nop}.txt')

        #     plt.title(f'Grad Free Energy for OP-{Nop}')
        #     im = plt.imshow( np.reshape(FE_grad[:,2], (Nv,Nu)), cmap='gist_heat' )
        #     plt.colorbar(im)
        #     plt.show()

    # elif Nu == 1: # basically for the 1D case
    #     # plot only the slices, since the 2D view is not useful here
    #     plt.title("Slices top to bottom")
    #     for i in range(Nop):
    #         slc = A[i][0]
    #         plt.plot( np.linspace(0,h*Nv,len(slc)), slc, label=f"A[{i}]" )
    #     plt.legend()
    #     plt.show()

# how the python interpreter will know
#   to run our function called 'main()'
if __name__ == "__main__":
    # call the main function with the arguments given in the terminal call
    main(sys.argv[1:])
