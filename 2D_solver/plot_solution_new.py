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
def plot_OP_comps_and_slices(Nop, organized_array, ext, h, size):

    # if it's going to be hard to see in the 3d plot...
    if Nop > 3:
        # plot the 2D solution
        for i in range(Nop):
            plt.xlabel(r'$z/\xi$')
            plt.ylabel(r'$x/\xi$')
            plt.title(f'OP component #{i}')
            im = plt.imshow(organized_array[i], extent=ext)
            plt.colorbar(im)
            plt.show()
            plt.clf()
        
    # plot 3D, all components
    plt.title(f'OP-{Nop}')
    ax2 = plt.subplot(111, projection='3d')
    ax2.set_xlabel(r'$x/\xi$')
    ax2.set_ylabel(r'$z/\xi$')
    ax2.set_zlabel(r'$|A_{\alpha i}|$')
    X, Y = np.meshgrid( np.linspace(ext[2],ext[3],len(organized_array[0][0])),
                        np.linspace(ext[0],ext[1],len(organized_array[0])) )
    for i in range(Nop):
        ax2.scatter(X,Y,organized_array[i],label=f'OP comp #{i+1}')
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
        exit()

    # read values from conditions file
    Nop, Nu, Nv, h, OP_array, labels = readSolutionFile(argv[0])

    # sort the OP_array for ease of plotting
    A = np.array([ np.reshape(OP_array[:,2*i+2], (Nv,Nu)) for i in range(Nop) ])

    # the domain extents for the imshow calls
    ext = [0*h, Nv*h, 0*h, Nu*h]

    # debugging: visualize the initial guess
    if debug: # TODO: check if this is working!
        # Initial guess -- for debugging
        plt.title("initial guess")
        # op = np.loadtxt(f'initGuess{Nop}.txt') # get the guess from the saved file from the C++ code in gl_fdm.cpp
        Nop, Nu, Nv, h, op, labels = readSolutionFile(f'initGuess{Nop}.txt')
        initOP = np.array([ np.reshape(op[:,2*i+2], (Nv,Nu)) for i in range(Nop) ])

        if Nu > 1:
            plot_OP_comps_and_slices(Nop, initOP, ext, h, Nv)
        elif Nu == 1: # basically for the 1D case
            # plot only the slices, since the 2D view is not useful here
            plt.title("Slices top to bottom")
            for i in range(Nop):
                slc = initOP[i][0]
                plt.plot( np.linspace(0,h*Nv,len(slc)), slc, label=f"initOP[{i}]" )
            plt.legend()
            plt.show()

    if Nu > 1:
        plot_OP_comps_and_slices(Nop, A, ext, h, Nv)

        if debug: # TODO: make sure this is working!
            # Plot the bulk free energy
            FE_bulk = np.loadtxt(f'bulkRHS_FE{Nop}.txt')

            plt.title(f'Bulk Free Energy for OP-{Nop}')
            im = plt.imshow( np.reshape(FE_bulk[:,2], (Nv,Nu)), extent=ext )
            plt.colorbar(im)
            plt.show()
            
            # Plot the grad free energy
            FE_grad = np.loadtxt(f'FEgrad{Nop}.txt')

            plt.title(f'Grad Free Energy for OP-{Nop}')
            im = plt.imshow( np.reshape(FE_grad[:,2], (Nv,Nu)) )
            plt.colorbar(im)
            plt.show()
        
        # but we'll always plot the total
        Ncols, Nu, Nv, h, totalFE, labels = read_FE_File(f'totalFE{Nop}.txt')
        # totalFE = np.loadtxt(f'totalFE{Nop}.txt')
        FE_on_grid = np.reshape(totalFE[:,2], (Nv,Nu))

        fig, (ax1,ax2) = plt.subplots(1,2)
        fig.suptitle(f'Total Free Energy for OP-{Nop}')

        ax1.set_xlabel(r'$z/\xi$')
        ax1.set_ylabel('FE')
        ax1.set_title('Total FE profile')
        ax1.plot(np.linspace(ext[0],ext[1],Nu), FE_on_grid[len(FE_on_grid)//2,:])

        ax2.set_xlabel(r'$z/\xi$')
        ax2.set_ylabel(r'$x/\xi$')
        im = ax2.imshow(FE_on_grid, extent=ext)
        # ax2.set_aspect(Nv/Nu)
        fig.colorbar(im,ax=ax2)
        plt.show()

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
    # readSolutionFile('solution5.txt')
