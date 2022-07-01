import sys # for arguments in terminal call
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def read_columns_from_file(file_name):
    column_dict = {} # holds all data from each column; the column label
                      #   is the key with value being an array of the data
    labels = []       # list for all labels
    Nu, Nv = -1, -1   # for grid size
    h = -1            # step size

    line_count = 0
    for line in open(file_name,'r'):
        if line_count == 0: # change to 1 if file has column numbers
            labels = line.split() # the first line should give us all the labels of the columns
            # start the lists in the dictionary
            for label in labels:
                column_dict[label] = []
        elif line[0] != '#':
            # get all the values from the line
            values = list(map(float,line.split()))

            if len(values) > 1:
                for i, label in enumerate(labels):
                    column_dict[label].append(values[i]) # put values into their columns
            else:
                # must be a blank line
                if Nv < 0: Nv = line_count-1 # we now know Nv!
                line_count -= 1 # decrement once since we'll count it later, but we only count lines with data on it
        line_count += 1
    # now that we know Nv and the total number of lines...
    Nu = line_count//Nv # MAY HAVE PROBLEMS FOR 1D SYSTEMS

    # calculate step size
    h = max( column_dict[labels[0]][1]-column_dict[labels[0]][0], column_dict[labels[1]][1]-column_dict[labels[1]][0] )
    
    # convert and reshape all items of the column dictionary
    if Nu > 1 and Nv > 1: # if it is a 2D system
        for label in labels:
            column_dict[label] = np.reshape(np.array(column_dict[label]), (Nu,Nv))
        
    return column_dict, labels # can also return Nu, Nv, h

# the main code to run
def main(argv): # argv will be like: [ file_name, [Nop] ]
    if len(argv) < 1:
        print("ERROR: you forgot to call this script with argumnets! I need a file_name.")
        exit()
    
    file_name = argv[0]
    Nop = 5 # UPDATE THIS!
    if len(argv) == 2: Nop == argv[1] # or get Nop from the command line
    CD, labels = read_columns_from_file(file_name) # can also get "Nu, Nv, h" from this function

    if input("Is this for a cylindrical system? (y/n): ") == 'y':
        custom_labels = [r'$A_{rr}$', r'$A_{\phi \phi}$', r'$A_{zz}$', r'$A_{zr}$', r'$A_{rz}$'] # cylindrical
        x_axis_label = rf'$r/\xi_0$'
    else:
        custom_labels = [r'$A_{xx}$', r'$A_{yy}$', r'$A_{zz}$', r'$A_{zx}$', r'$A_{xz}$']
        x_axis_label = rf'$x/\xi_0$'
    z_axis_label = rf'$z/\xi_0$'

    clr_map = LinearSegmentedColormap.from_list("", ["navy","blue","white","red","maroon"]) # colors for the gradient in plots; low value to high
    def PColorMeshPlot(v_col, lower_lim, upper_lim, ax, cmap=None, uShift=0, vShift=0, u_col=0, w_col=1):
        if cmap == None:
            cmap = ListedColormap(clr_map(np.linspace(
                # play with these values! must be between 0.0 and 1.0, inclusive
                lower_lim, # the lower limit/percent of the color gradient
                upper_lim, # the upper limit/percent
                128  # don't change! must be 128
            )))
        else:
            cmap = ListedColormap(cm.get_cmap(cmap,128)(np.linspace(
                # play with these values! must be between 0.0 and 1.0, inclusive
                lower_lim, # the lower limit/percent of the color gradient
                upper_lim, # the upper limit/percent
                128  # don't change! must be 128
            )))
        return ax.pcolormesh(
            CD[labels[u_col]]-uShift, # u-data, from column u_col # this one shouldn't need to be changed
            CD[labels[w_col]]-vShift, # w-data, from column w_col # this one shouldn't need to be changed
            CD[labels[v_col]], # v-data, from column v_col
            shading='gouraud', # to make the plot look more smooth
            cmap=cmap
        )

    OP_axs = [None]*Nop
    FE_ax, grad_FE_ax, FE_ref_ax, empty_ax2 = None, None, None, None

    # shift values...for x- and y-axes
    u_shift = 0
    v_shift = 0


    # SET UP PLOTTING FIGURES
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # for 3-component OP plotting
    if Nop == 3:
        fig, axes = plt.subplots(2,3)
        ((OP_axs[0], OP_axs[1], OP_axs[2]),
        (grad_FE_ax, FE_ax,     FE_ref_ax)) = axes # unpack the axes

    # for 5-component OP plotting
    if Nop == 5:
        fig, axes = plt.subplots(3,3)
        ((grad_FE_ax, OP_axs[0], OP_axs[3]),
        (FE_ax,       OP_axs[1], OP_axs[4]),
        (FE_ref_ax,   OP_axs[2], empty_ax2)) = axes # unpack the axes
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # SET ORDER PARAMETER TITLES FOR SUBPLOTS
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # for both 3 & 5
    for i, ax in enumerate(OP_axs): # set labels for the OP components
        ax.set_title(custom_labels[i])
    # turn off unused axes
    # FE_ref_ax.axis('off') # BUT DON'T TURN THIS ONE OFF IF YOU WANT TO PLOT SOMETHING ELSE...LIKE A REFERENCE ENERGY DENSITY!
    if empty_ax2 != None: empty_ax2.axis('off')
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # OP COMPONENT 1
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    if Nop == 5: OP_axs[0].axes.yaxis.set_ticks([]) # for 5 comp
    if Nop == 3: OP_axs[0].set_ylabel(z_axis_label) # for 3 comp
    OP_axs[0].axes.xaxis.set_ticks([])
    pcm = PColorMeshPlot(2, 0.7, 1.0, OP_axs[0], uShift=u_shift, vShift=v_shift) # change these values!
    plt.colorbar(pcm,ax=OP_axs[0]) # show the colorbar for this 2D plot
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # OP COMPONENT 2
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    OP_axs[1].axes.yaxis.set_ticks([])
    OP_axs[1].axes.xaxis.set_ticks([])
    pcm = PColorMeshPlot(4, 0.7, 1.0, OP_axs[1], uShift=u_shift, vShift=v_shift) # change these values!
    plt.colorbar(pcm,ax=OP_axs[1])
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # OP COMPONENT 3
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    OP_axs[2].set_xlabel(x_axis_label)
    # if Nop == 5: OP_axs[2].set_ylabel(z_axis_label) # for 5 comp
    OP_axs[2].axes.yaxis.set_ticks([])
    if Nop == 3: OP_axs[2].axes.yaxis.set_ticks([]) # for 3 comp
    pcm = PColorMeshPlot(6, 0.25, 0.8, OP_axs[2], uShift=u_shift, vShift=v_shift) # change these values!
    plt.colorbar(pcm,ax=OP_axs[2])
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # OP COMPONENT 4
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    if Nop == 5:
        OP_axs[3].axes.xaxis.set_ticks([])
        OP_axs[3].axes.yaxis.set_ticks([])
        pcm = PColorMeshPlot(8, 0.4, 0.7, OP_axs[3], uShift=u_shift, vShift=v_shift) # change these values!
        plt.colorbar(pcm,ax=OP_axs[3])
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    

    # OP COMPONENT 5
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    if Nop == 5:
        OP_axs[4].axes.yaxis.set_ticks([])
        OP_axs[4].set_xlabel(x_axis_label)
        pcm = PColorMeshPlot(10, 0.4, 0.6, OP_axs[4], uShift=u_shift, vShift=v_shift) # change these values!
        plt.colorbar(pcm,ax=OP_axs[4])
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    

    # GRADIENT FREE ENERGY
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    grad_FE_ax.set_title(r'$FE_{Grad}$')
    grad_FE_ax.set_ylabel(z_axis_label)
    if Nop == 5: grad_FE_ax.axes.xaxis.set_ticks([]) # for 5 comp
    if Nop == 3: grad_FE_ax.set_xlabel(x_axis_label) # for 3 comp
    # change these values! # the first value here will have to be different for 3 & 5 comp OP: 12-15 for 5comp, 8-11 for 3comp
    pcm = PColorMeshPlot(14, 0.5, 1., grad_FE_ax, uShift=u_shift, vShift=v_shift, cmap='gist_heat')
    plt.colorbar(pcm,ax=grad_FE_ax)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # TOTAL FREE ENERGY
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    FE_ax.set_title(r'$FE_{Total}$')
    if Nop == 5: FE_ax.set_ylabel(z_axis_label) # for 5 comp
    # FE_ax.set_xlabel(x_axis_label)
    FE_ax.axes.xaxis.set_ticks([])
    if Nop == 3: FE_ax.axes.yaxis.set_ticks([]) # for 3 comp
    # change these values! # the first value here will have to be different for 3 & 5 comp OP: 12-15 for 5comp, 8-11 for 3comp
    pcm = PColorMeshPlot(12, 0., 0.5, FE_ax, uShift=u_shift, vShift=v_shift, cmap='gist_heat')
    plt.colorbar(pcm,ax=FE_ax)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # REFERENCE ENERGY
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    FE_ref_ax.set_title(r'$FE - FE_B$')
    if Nop == 5: FE_ref_ax.set_ylabel(z_axis_label) # for 5 comp
    FE_ref_ax.set_xlabel(x_axis_label)
    if Nop == 3: FE_ref_ax.axes.yaxis.set_ticks([]) # for 3 comp
    # change these values! # the first value here will have to be different for 3 & 5 comp OP: 12-15 for 5comp, 8-11 for 3comp
    pcm = PColorMeshPlot(15, 0.2, 1., FE_ref_ax, uShift=u_shift, vShift=v_shift, cmap='PuOr_r')
    plt.colorbar(pcm,ax=FE_ref_ax)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # THE BASIC WAY TO JUST PLOT ONE FIGURE
    # plt.clf() # make sure the figure is clear
    # plt.title('Title')
    # plt.xlabel(x_axis_label)
    # plt.ylabel(z_axis_label)
    # pcm = PColorMeshPlot(2,0.5,0.9,plt)
    # plt.colorbar(pcm) # show the colorbar for this 2D plot
    # plt.show()
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


    plt.show() # show the finished figure

# how the python interpreter will know to run our function called 'main()'
if __name__ == "__main__":
    main(sys.argv[1:]) # call the main function with the arguments given in the terminal call
