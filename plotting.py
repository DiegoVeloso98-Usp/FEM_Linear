import numpy as np
import matplotlib.pyplot as plt
import scienceplots
import solver
from matplotlib.tri import Triangulation

plt.style.use(['science','notebook','grid'])
plt.figure(figsize=(15, 6))  # Set the size of the plot area

def undeformedPlotting(ELEMS,NODES):


    for elem in ELEMS:
        # Remove last 5 columns
        elem = elem[:-5]
        
        x_coords = np.array([])
        y_coords = np.array([])

        for node in elem:
            # Node indices for this element 
            node_indices = int(node - 1)  
            
            # Coordinates of these nodes 
        
            x_coords=np.append(x_coords,NODES[node_indices, 0])
            y_coords=np.append(y_coords,NODES[node_indices, 1])
        
        # Closing the loop by adding the first node at the end

        x_coords=np.append(x_coords,x_coords[0])
        y_coords=np.append(y_coords,y_coords[0])
        

        plt.fill(x_coords, y_coords, 'lightblue', edgecolor='black', alpha=0.5)
        plt.plot(x_coords, y_coords, color='black', linewidth=0.5)
    return

def amplifiedDeformed(NODES,nnodes, u,a):
    displa=np.copy(NODES)
    for in_node in range(nnodes):
        displa[in_node, 0] += a * u[2 * in_node]
        displa[in_node, 1] += a * u[2 * in_node + 1]
    return displa

def deformedPlotting(ELEMS,NODES, u, nnodes,a):

    for elem in ELEMS:
        elem = elem[:-5]

        displa=amplifiedDeformed(NODES,nnodes, u,a)

        x_coords = np.array([])
        y_coords = np.array([])

        for node in elem:
            node_indices = int(node - 1)  
            
            x_coords=np.append(x_coords,displa[node_indices, 0]) 
            y_coords=np.append(y_coords,displa[node_indices, 1])
        
        x_coords=np.append(x_coords,x_coords[0]) 
        y_coords=np.append(y_coords,y_coords[0])

        
        plt.fill(x_coords, y_coords, 'lightgreen', edgecolor='black', alpha=0.5)


    
    legend1 = plt.legend(loc='upper right', fontsize=1,bbox_to_anchor=(.9, .95),title="Deformed")
    legend1.get_frame().set_facecolor('lightgreen')

    legend2 = plt.legend(loc='upper right',fontsize=1, bbox_to_anchor=(.9, .87), title="Undeformed")
    legend2.get_frame().set_facecolor('lightblue')
    
    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)

    plt.xlim(min(NODES[:, 0])-.5 , max(NODES[:, 0]) + .5)
    plt.ylim(min(NODES[:, 1]) -.5, max(NODES[:, 1])+.5 )
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Filled Elements Plot')
    plt.grid(True)
 
    plt.show()


    return
'''
def deformedPlotting(ELEMS, NODES, SIGMA, u, nnodes, a):
    displa = amplifiedDeformed(NODES, nnodes, u, a)  # Calculate the deformed positions

    for elem in ELEMS:
        # Correctly remove last 5 columns
        elem = elem[:-5]

        # Initialize lists to store coordinates and corresponding SIGMA values for the current element
        x_coords = []
        y_coords = []
        sigma_values = []

        for node in elem:
            # Extract node indices for this element (convert from 1-based to 0-based index)
            node_indices = int(node - 1)

            # Get the deformed coordinates of these nodes
            x_coords.append(displa[node_indices, 0])
            y_coords.append(displa[node_indices, 1])

            # Get the corresponding SIGMA value (use the first column of SIGMA)
            sigma_values.append(SIGMA[node_indices, 0])

        # Close the loop by adding the first node at the end
        # x_coords.append(x_coords[0])
        # y_coords.append(y_coords[0])

        # Plot the triangular element with a gradient based on SIGMA values
        triang = Triangulation(x_coords, y_coords)
        plt.tripcolor(triang, sigma_values, shading='gouraud', cmap='viridis')

    # Add legends and labels
    legend1 = plt.legend(loc='upper right', fontsize=1, bbox_to_anchor=(.9, .95), title="Deformed")
    legend1.get_frame().set_facecolor('lightgreen')

    legend2 = plt.legend(loc='upper right', fontsize=1, bbox_to_anchor=(.9, .87), title="Undeformed")
    legend2.get_frame().set_facecolor('lightblue')

    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)

    # Set plot limits and labels
    plt.xlim(min(NODES[:, 0]) - 0.5, max(NODES[:, 0]) + 0.5)
    plt.ylim(min(NODES[:, 1]) - 0.5, max(NODES[:, 1]) + 0.5)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Gradient Plot of Deformed Elements')
    plt.grid(True)
    
    # Add color bar to reflect the gradient values
    plt.colorbar(label='Sigma Values')

    plt.show()

    return

'''
def displacedParticles(NODESP, U_P,NODES,U,a):

    # Calculate the displaced positions
    U = np.array(U)
    U_P = np.array(U_P)
    disp_part = NODESP + a * U_P.reshape(-1, 2)
    disp_nodes = NODES + a * U.reshape(-1, 2)

    # Create the plot
    plt.figure(figsize=(8, 8))

    plt.scatter(NODESP[:, 0], NODESP[:, 1], color='blue', label='Initial Particle Nodes', s=3)

    plt.scatter(disp_part[:, 0], disp_part[:, 1], color='red', label='Displaced Particle Nodes', s=3)

    plt.scatter(NODES[:, 0], NODES[:, 1], color='black', label='Structure Nodes', s=3)
    plt.scatter(disp_nodes[:, 0], disp_nodes[:, 1], color='green', label='Displaced Structure Nodes', s=3)
    
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title('Initial and Displaced Node Positions')
    plt.legend()
    plt.grid(True)
    
    plt.axis('equal')
    
    plt.show()
    return
