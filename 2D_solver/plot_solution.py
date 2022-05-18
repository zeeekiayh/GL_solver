import sys # for arguments in terminal call
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# get parameters from f'conditions{Nop}.txt'
# read in the conditions from the file
def readConditions(Nop):
    conditions = []
    line_count = 0
    for line in open(f'conditions{Nop}.txt','r'):
        num = list(map(str, line.split()))  # put all the values into a list
        if line_count == 0:
            conditions.append(float(num[0]))
        elif line_count == 13+2*(Nop-3): # any better way to do this?! to account for more lines in the conditions file for large OP?
            conditions.append(float(num[0]))
            conditions.append(float(num[1]))
        elif line_count == 14+2*(Nop-3): # any better way to do this?! to account for more lines in the conditions file for large OP?
            conditions.append(float(num[0]))
        line_count+=1

    # print(f'{conditions = }')
    # define variables from the conditions list
    opSize = round(conditions[0])
    size_x = round(conditions[1])
    size_z = round(conditions[2])
    step   = conditions[3]
    if not opSize == Nop: print(f"WARNING: this python3 code has not been written to handle {opSize} components.")

    return opSize, size_x, size_z, step

# plot all the OP components in 2D and slices
def plot_OP_comps_and_slices(Nop, organized_array, ext, h, size):
    # plot the 2D solution
    for i in range(Nop):
        plt.title(f'OP component #{i}')
        im = plt.imshow(organized_array[i], extent=ext)
        plt.colorbar(im)
        plt.show()
        plt.clf()
        
    # plot 3D, all components
    plt.title(f'OP component #{i}')
    ax2 = plt.subplot(111, projection='3d')
    X, Y = np.meshgrid( np.linspace(ext[0],ext[1],len(organized_array[i][0])),
                        np.linspace(ext[2],ext[3],len(organized_array[i])) )
    surf = None
    for i in range(Nop):
        surf = ax2.plot_surface(X,Y,organized_array[i])
    plt.colorbar(surf)
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

    fig = plt.figure(figsize=(20,20))
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
    ext = [min(OP_array[:,0]), max(OP_array[:,0]), min(OP_array[:,1]), max(OP_array[:,1])]

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
    main(sys.argv[1:])
