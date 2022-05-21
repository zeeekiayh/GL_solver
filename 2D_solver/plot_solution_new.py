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
            print(f'{labels = }')
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

    Nop = (len(labels)-2)/2 # take off the 2 cols for location; /2 b/c of complex cols

    # print(f'{u_0_count = }')
    # print(f'{v_0_count = }')
    Nu = 1+v_0_count # YES! the u and v are switched! they need to be!
    Nv = 1+u_0_count
        
    return Nop, Nu, Nv, h, np_data_array, labels

def read_FE_File(file_name):
    np_data_array = [] # to be converted to np.array() later
    labels = [] # to store the labels of all the columns
    Ncols = 0
    h = 0 # step size; we'll assume it will be the same in all directions

    line_count = 0
    for line in open(file_name,'r'):
        # the first line should give us all the labels of the columns
        if line_count == 0:
            labels = line.split()
            Ncols = len(labels)
            np_data_array.append([]*Ncols)
            print(f'{labels = }')
            print(f'{np_data_array = }')
        else:
            # start reading in the values into np_data_array
            pass
        line_count += 1
        
    # return opSize, size_x, size_z, step
    return Nop, Nu, Nv, h, np_data_array, labels

# plot all the OP components in 2D and slices
def plot_OP_comps_and_slices(Nop, organized_array, ext, h, size):
    plt.xlabel(r'$z/\xi$')
    plt.ylabel(r'$x/\xi$')

    # plot the 2D solution
    for i in range(Nop):
        plt.title(f'OP component #{i}')
        im = plt.imshow(organized_array[i], extent=ext)
        plt.colorbar(im)
        plt.show()
        plt.clf()
        
    # plot 3D, all components
    plt.title(f'OP-{Nop}')
    ax2 = plt.subplot(111, projection='3d')
    ax2.set_xlabel(r'$z/\xi$')
    ax2.set_ylabel(r'$x/\xi$')
    ax2.set_zlabel(r'$|A_{\alpha i}|$')
    X, Y = np.meshgrid( np.linspace(ext[2],ext[3],len(organized_array[0][0])),
                        np.linspace(ext[0],ext[1],len(organized_array[0])) )
    surf = None
    for i in range(Nop):
        surf = ax2.plot_surface(X,Y,organized_array[i],label=f'OP comp #{i+1}')
        surf._facecolors2d = surf._facecolor3d # We need these because of a strange bug in matplotlib;
        surf._edgecolors2d = surf._edgecolor3d # see https://stackoverflow.com/questions/55531760/is-there-a-way-to-label-multiple-3d-surfaces-in-matplotlib
    plt.legend()
    plt.show()
    plt.clf()

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

    # fig = plt.figure(figsize=(20,20))
    # gs = gridspec.GridSpec()

    # read in Nop from arguments
    Nop = int(argv[0])

    # read values from conditions file
    opSize, size_x, size_z, step = readConditions(Nop)

    # load the OP values from the solution file
    OP_array = np.loadtxt(f'solution{Nop}.txt')

    # sort the OP_array for ease of plotting
    A = np.array([ np.reshape(OP_array[:,2*i+2], (size_z,size_x)) for i in range(Nop) ])

    # the domain extents for the imshow calls
    ext = [min(OP_array[:,1]), max(OP_array[:,1]), min(OP_array[:,0]), max(OP_array[:,0])]

    # debugging: visualize the initial guess
    if debug:
        # Initial guess -- for debugging
        plt.title("initial guess")
        op = np.loadtxt(f'initGuess{Nop}.txt') # get the guess from the saved file from the C++ code in gl_fdm.cpp
        initOP = np.array([ np.reshape(op[:,2*i+2], (size_z,size_x)) for i in range(Nop) ])

        if size_x > 1:
            plot_OP_comps_and_slices(Nop, initOP, ext, step, size_z)
        elif size_x == 1: # basically for the 1D case
            # plot only the slices, since the 2D view is not useful here
            plt.title("Slices top to bottom")
            for i in range(Nop):
                slc = initOP[i][0]
                plt.plot( np.linspace(0,step*size_z,len(slc)), slc, label=f"initOP[{i}]" )
            plt.legend()
            plt.show()

    if size_x > 1:
        plot_OP_comps_and_slices(Nop, A, ext, step, size_z)

        if debug:
            # Plot the bulk free energy
            FE_bulk = np.loadtxt(f'bulkRHS_FE{Nop}.txt')

            plt.title(f'Bulk Free Energy for OP-{Nop}')
            im = plt.imshow(np.reshape(FE_bulk[:,2], (size_z,size_x)), extent=ext)
            plt.colorbar(im)
            plt.show()
            
            # Plot the grad free energy # CONTINUE HERE!
            FE_grad = np.loadtxt(f'FEgrad{Nop}.txt')

            plt.title(f'Grad Free Energy for OP-{Nop}')
            im = plt.imshow(np.reshape(FE_grad[:,2], (size_z,size_x)))
            plt.colorbar(im)
            plt.show()
        
        # but we'll always plot the total
        FE_bulk = np.loadtxt(f'totalFE{Nop}.txt')
        FE_on_grid = np.reshape(FE_bulk[:,2], (size_z,size_x))

        fig, (ax1,ax2) = plt.subplots(1,2)
        fig.suptitle(f'Total Free Energy for OP-{Nop}')

        ax1.set_ylabel(r'$z/\xi$')
        ax1.set_xlabel(r'$x/\xi$')
        ax1.set_title('Total FE profile')
        ax1.plot(np.linspace(ext[0],ext[1],size_x), FE_on_grid[len(FE_on_grid)//2,:])

        ax2.set_xlabel(r'$z/\xi$')
        ax2.set_ylabel(r'$x/\xi$')
        im = ax2.imshow(FE_on_grid, extent=ext)
        ax2.set_aspect(size_z/size_x)
        fig.colorbar(im,ax=ax2)
        plt.show()

    elif size_x == 1: # basically for the 1D case
        # plot only the slices, since the 2D view is not useful here
        plt.title("Slices top to bottom")
        for i in range(Nop):
            slc = A[i][0]
            plt.plot( np.linspace(0,step*size_z,len(slc)), slc, label=f"A[{i}]" )
        plt.legend()
        plt.show()

# how the python interpreter will know
#   to run our function called 'main()'
if __name__ == "__main__":
    # call the main function with the arguments given in the terminal call
    # main(sys.argv[1:])
    readSolutionFile('solution5.txt')
