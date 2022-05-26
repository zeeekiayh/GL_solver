import sys # for arguments in terminal call
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
            # print(f'{labels = }')
        else:
            # start reading in the values into np_data_array
            np_data_array.append(list(map(float,line.split())))

            if (abs(np_data_array[-1][0]) < 1e-8):
                u_0_count += 1
            if (abs(np_data_array[-1][1]) < 1e-8):
                v_0_count += 1
        line_count += 1
    # print(f'{np_data_array = }')

    # convert the data array
    np_data_array = np.array(np_data_array)
    # this is now in the exact for that we had in the other python plotting file!

    # step size; we'll assume it will be the same in all directions
    h = np_data_array[1,0] - np_data_array[0,0]
    if abs(h) < 1e-8: h = np_data_array[1,1] - np_data_array[0,1] # in case position values are in a different order
    # print(f'{h = }')

    Nop = (len(labels)-2)//2 # take off the 2 cols for location; /2 b/c of complex cols

    # print(f'{u_0_count = }')
    # print(f'{v_0_count = }')
    Nu = int(v_0_count) # YES! the u and v are switched! they need to be!
    Nv = int(u_0_count)
        
    return Nop, Nu, Nv, h, np_data_array, labels

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

# plot all the OP components in 2D and slices
def plot_OP_comps_and_slices(Nop, organized_array, X, Z, ext, labels, FE_file_name):
    # initialize all axes
    axs = []
    for _ in range(Nop): axs.append(None)
    # other axes we may use
    FE_ax, FE_prof_ax, ax_empty, ax3D_if_wanted, fig = None, None, None, None, None

    # shape the plot based on OP size
    if Nop == 3:   fig, ((axs[0], axs[1], axs[2]), (FE_prof_ax, FE_ax, ax3D_if_wanted)) = plt.subplots(3,2)
    elif Nop == 5: fig, ((axs[0], axs[1], axs[2]), (axs[3], axs[4], ax_empty), (FE_prof_ax, FE_ax, ax3D_if_wanted)) = plt.subplots(3,3)
    else: print(f"Implement 'plot_OP_comps_and_slices' for {Nop = }.")

    fig.suptitle(f'OP-{Nop}')

    # plot the 2D solution
    for i in range(Nop):
        if i == 0:
            axs[i].set_xlabel(r'$z/\xi_0$', labelpad=-15)
            axs[i].set_ylabel(r'$x/\xi_0$', labelpad=0)
            axs[i].set_title(f'comp: {labels[i]}', y=1.0)
        elif i > 2:
            axs[i].axes.xaxis.set_ticks([])
            axs[i].axes.yaxis.set_ticks([])
            axs[i].set_ylabel(f'comp: {labels[i]}')
        elif i > 0:
            axs[i].set_title(f'comp: {labels[i]}', y=1.0)
            axs[i].axes.xaxis.set_ticks([])
            axs[i].axes.yaxis.set_ticks([])
        
        # don't show the blank plots
        if ax_empty != None: ax_empty.axis('off')
        if ax3D_if_wanted != None: ax3D_if_wanted.axis('off')

        im = axs[i].imshow(organized_array[i], extent=ext, cmap='bwr')
        plt.colorbar(im,ax=axs[i])

    # get the FE data
    Ncols, Nu, Nv, h, totalFE, labels = read_FE_File(FE_file_name)
    FE_on_grid = np.reshape(totalFE[:,2], (Nv,Nu))

    # plot the FE profile
    FE_prof_ax.set_xlabel(r'$z/\xi_0$', labelpad=0)
    FE_prof_ax.set_ylabel('FE', labelpad=0)
    FE_prof_ax.set_title('Total FE profile', y=1.0)
    FE_prof_ax.plot(np.linspace(ext[2],ext[3],Nu), FE_on_grid[len(FE_on_grid)//2,:])

    # plot the 2D heatmap of the FE
    FE_ax.set_xlabel(r'$z/\xi_0$', labelpad=0)
    FE_ax.set_ylabel(r'$x/\xi_0$', labelpad=0)
    FE_ax.set_title('Total FE', y=1.0)
    im = FE_ax.imshow(FE_on_grid.transpose(), extent=[0,ext[1],0,ext[3]], cmap='gist_heat')
    fig.colorbar(im,ax=FE_ax)

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
    #     # ax.scatter(X,Z,organized_array[i],label=f'OP comp {labels[i]}')
    #     # If the 3D plot is too crowded, use thinned out arrays!
    #     ax.scatter(X[::2,::2],Z[::2,::2],organized_array[i][::2,::2],label=f'OP comp #{labels[i]}')
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

    # read values from conditions file
    Nop, Nu, Nv, h, OP_array, labels = readSolutionFile(argv[0])

    # sort the OP_array for ease of plotting
    A = np.array([ np.reshape(OP_array[:,2*i+2], (Nv,Nu)).transpose() for i in range(Nop) ])
    X, Z = np.reshape(OP_array[:,0], (Nv,Nu)), np.reshape(OP_array[:,1], (Nv,Nu))
    X = X.transpose() # make sure they are good for display!
    Z = Z.transpose()

    # the domain extents for the imshow calls
    ext = [0*h, Nv*h, 0*h, Nu*h]

    # debugging: visualize the initial guess
    if debug:
        # Initial guess -- for debugging
        plt.title("initial guess")
        # get the guess from the saved file from the C++ code in gl_fdm.cpp
        Nop, Nu, Nv, h, op, labels = readSolutionFile(f'initGuess{Nop}.txt')
        initOP = np.array([ np.reshape(op[:,2*i+2], (Nv,Nu)) for i in range(Nop) ])

        if Nu > 1:
            plot_OP_comps_and_slices( Nop, initOP, X, Z, ext, list(map(lambda l: l+"_guess", labels[2::2])), None )
        elif Nu == 1: # basically for the 1D case
            # plot only the slices, since the 2D view is not useful here
            plt.title("Slices top to bottom")
            for i in range(Nop):
                slc = initOP[i][0]
                plt.plot( np.linspace(0,h*Nv,len(slc)), slc, label=f"initOP[{i}]" )
            plt.legend()
            plt.show()

    if Nu > 1:
        plot_OP_comps_and_slices(Nop, A, X, Z, ext, labels[2::2], f'totalFE{Nop}.txt' )

        if debug:
            # Plot the bulk free energy
            Ncols, Nu, Nv, h, FE_bulk, labels = read_FE_File(f'bulkRHS_FE{Nop}.txt')

            plt.title(f'Bulk Free Energy for OP-{Nop}')
            im = plt.imshow( np.reshape(FE_bulk[:,2], (Nv,Nu)), extent=ext, cmap='gist_heat' )
            plt.colorbar(im)
            plt.show()
            
            # Plot the grad free energy
            Ncols, Nu, Nv, h, FE_grad, labels = read_FE_File(f'FEgrad{Nop}.txt')

            plt.title(f'Grad Free Energy for OP-{Nop}')
            im = plt.imshow( np.reshape(FE_grad[:,2], (Nv,Nu)), cmap='gist_heat' )
            plt.colorbar(im)
            plt.show()
        
        # # but we'll always plot the total
        # Ncols, Nu, Nv, h, totalFE, labels = read_FE_File(f'totalFE{Nop}.txt')
        # FE_on_grid = np.reshape(totalFE[:,2], (Nv,Nu))

        # fig, (ax1,ax2) = plt.subplots(2,1)
        # fig.suptitle(f'Total Free Energy for OP-{Nop}')

        # ax1.set_xlabel(r'$z/\xi_0$')
        # ax1.set_ylabel('FE')
        # ax1.set_title('Total FE profile')
        # ax1.plot(np.linspace(ext[2],ext[3],Nu), FE_on_grid[len(FE_on_grid)//2,:])

        # ax2.set_xlabel(r'$x/\xi_0$')
        # ax2.set_ylabel(r'$z/\xi_0$')
        # im = ax2.imshow(FE_on_grid, extent=[0,ext[3],0,ext[1]], cmap='gist_heat')
        # # ax2.set_aspect(Nv/Nu)
        # fig.colorbar(im,ax=ax2)
        # plt.show()

    elif Nu == 1: # basically for the 1D case
        # plot only the slices, since the 2D view is not useful here
        plt.title("Slices top to bottom")
        for i in range(Nop):
            slc = A[i][0]
            plt.plot( np.linspace(0,h*Nv,len(slc)), slc, label=f"A[{i}]" )
        plt.legend()
        plt.show()

# how the python interpreter will know
#   to run our function called 'main()'
if __name__ == "__main__":
    # call the main function with the arguments given in the terminal call
    main(sys.argv[1:])