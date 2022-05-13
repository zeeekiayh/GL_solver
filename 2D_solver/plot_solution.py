import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

################################
# Initial guess -- for debugging
# plt.title("initial guess")

# op = np.loadtxt(f'initGuess{Nop}.txt')

# Axx = np.reshape(op[:,2], (size_x,size_z))
# im = plt.imshow(Axx)
# plt.colorbar(im)
# plt.show()

# Ayy = np.reshape(op[:,4], (size_x,size_z))
# im = plt.imshow(Ayy)
# plt.colorbar(im)
# plt.show()

# Azz = np.reshape(op[:,6], (size_x,size_z))
# im = plt.imshow(Azz)
# plt.colorbar(im)
# plt.show()

# plt.clf()

# plt.title("Guess: slices top to bottom")
# plt.plot(Axx[:,len(Axx)//2],label="Axx")
# plt.plot(Ayy[:,len(Ayy)//2],label="Ayy")
# plt.plot(Azz[:,len(Azz)//2],label="Azz")
# plt.legend()
# plt.show()

def readConditions(file_name):
    # get parameters from f'conditions{Nop}.txt'
    # read in the conditions from the file
    conditions = []
    line_count = 0
    for line in open(file_name):
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

def main(argv):
    # argv should be like:
    #   [ Nop, [debug=false] ]
    if (len(argv)==2): print("Using debugging in python code...")
    elif (len(argv)!=1): print("Incorrect number of arguments in python call!")

    fig = plt.figure(figsize=(20,20))
    # gs = gridspec.GridSpec()

    # read in Nop from arguments
    Nop = int(argv[0])

    # read values from conditions file
    opSize, size_x, size_z, step = readConditions(f'conditions{Nop}.txt')

    # load the OP values from the solution file
    OP_array = np.loadtxt(f'solution{Nop}.txt')

    A = np.array([ np.reshape(OP_array[:,2*i+2], (size_x,size_z)) for i in range(Nop) ])
    # print(f'{A = }')
    if size_x > 1:
        for i in range(Nop):
            plt.title(f'OP component #{i}')
            im = plt.imshow(A[i], extent=[min(OP_array[:,0]),max(OP_array[:,0]), min(OP_array[:,1]), max(OP_array[:,1])])
            plt.colorbar(im)
            plt.show()

        # plot slices
        plt.clf()
        plt.title("Slices top to bottom")
        for i in range(Nop):
            slc = A[i][:,len(A[i])//2]
            plt.plot( np.linspace(0,step*size_z,len(slc)), slc, label=f"A[{i}]" )
        plt.legend()
        plt.show()
    elif size_x == 1: # basically for the 1D case
        plt.title("Slices top to bottom")
        for i in range(Nop):
            slc = A[i][0]
            plt.plot( np.linspace(0,step*size_z,len(slc)), slc, label=f"A[{i}]" )
        plt.legend()
        plt.show()

    #############################
    # Plot the bulk free energy #
    #############################
    FE_bulk = np.loadtxt(f'bulkRHS_FE{Nop}.txt')
    # print(f'{FE_bulk = }')
    for i in range(Nop):
        plt.title(f'Free Energy for OP-{Nop}; comp#{i}')
        im = plt.imshow(np.reshape(FE_bulk[:,2+2*i], (size_x,size_z)), extent=[min(OP_array[:,0]),max(OP_array[:,0]), min(OP_array[:,1]), max(OP_array[:,1])])
        plt.colorbar(im)
        plt.show()

    #############################
    # Plot the grad free energy # CONTINUE HERE!
    #############################
    # FE_grad = np.loadtxt(f'gradFE{Nop}.txt')
    # # print(f'{FE_grad = }')
    # for i in range(Nop):
    #     plt.title(f'Free Energy for OP-{Nop}; comp#{i}')
    #     im = plt.imshow(np.reshape(FE_grad[:,2+2*i], (size_x,size_z)))
    #     plt.colorbar(im)
    #     plt.show()

if __name__ == "__main__":
    # call the main function with the arguments given in the terminal call
    main(sys.argv[1:])